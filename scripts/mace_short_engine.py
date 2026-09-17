"""Exact pinned local MACE readout and its own Cartesian derivative."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_hybrid import EV_TO_KCAL,accepted_attempt,check_atoms
from mace_global_benchmark import snapshot,numerical_parent_gate
from mace_curvature import actual_collection,dft_source

SCHEMA='alquemia.mace_short_engine.v1'
COMPONENT='interaction_energy'
ADAPTER='polar_scale_shift_exact_readout_v1'
TOL={'energy_eV':1e-6,'energy_kcal_mol':.01,'force_max_eV_A':.001,
     'odd_absolute_kcal_mol':.02,'odd_relative':.05}


def evaluate(calc,atoms):
    """Stop the unmodified model after scale_shift, before global charge updates.

    Retain the original autograd graph and differentiate the captured scalar.
    No total-force reuse, independent local reimplementation or density output.
    """
    import torch
    if len(calc.models)!=1 or calc.use_compile or any(atoms.pbc):
        raise InvalidArtifact('single uncompiled nonperiodic model required')
    model=calc.models[0]
    if not hasattr(model,'lr_source_maps') or not hasattr(model,'scale_shift'):
        raise InvalidArtifact('pinned PolarMACE local readout unavailable')
    batch=calc._atoms_to_batch(atoms);data=batch.to_dict()
    if data['positions'].shape!=(len(atoms),3) or int(data['ptr'].numel())!=2:
        raise InvalidArtifact('padded/multiple-graph short evaluation unsupported')

    class ReadoutReady(Exception):
        def __init__(self,tensor):self.tensor=tensor

    def capture(_module,_inputs,output):
        raise ReadoutReady(output)

    hook=model.scale_shift.register_forward_hook(capture)
    try:
        try:
            model(data,compute_force=False,compute_stress=False,training=False)
        except ReadoutReady as caught:
            node_values=caught.tensor
        else:
            raise InvalidArtifact('model did not reach declared early-exit readout')
    finally:
        hook.remove()
    if node_values.shape!=(len(atoms),) or not data['positions'].requires_grad:
        raise InvalidArtifact('short readout shape/coordinate graph mismatch')
    value=node_values.sum()*calc.energy_units_to_eV
    gradient=torch.autograd.grad(value,data['positions'],create_graph=False,retain_graph=False)[0]
    return float(value.detach().cpu()),(-gradient/calc.length_units_to_A).detach().cpu().numpy()


def prepare(full,core,agreement,output):
    fc,fm=actual_collection(full);cc,cm=actual_collection(core)
    if (fm['checkpoint_label']!='medium' or cm['checkpoint_label']!='medium' or
            fm['software']!=cm['software'] or fm['model']['checkpoint']!=cm['model']['checkpoint'] or
            len(fm['tasks'])!=14 or len(cm['tasks'])!=20):
        raise InvalidArtifact('exact completed medium full and GGR curvature panels required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    pins=snapshot(out,(*fm['implementation'],*cm['implementation'],'mace_short_engine.py'))
    model=copy.deepcopy(fm['model']);model.update(energy_component=COMPONENT,short_adapter=ADAPTER,
                                               preparation_policy='exact_saved_full_and_GGR_curvature_v1')
    tasks=[]
    for prefix,c,m,path in (('full',fc,fm,full),('core',cc,cm,core)):
        for old in m['tasks']:
            t=copy.deepcopy(old);t.pop('cache_key');t['task_id']=prefix+'__'+old['task_id']
            t.update(source_task_id=old['task_id'],source_collection=record(path),energy_component=COMPONENT,
                     reference_component_energy_eV=c['rows'][old['task_id']]['energy_components_eV'][COMPONENT])
            tasks.append(t)
    m={'schema_version':SCHEMA,'protocol_id':'mace_polar_1m_short_energy_gradient_v1','agreement':record(agreement),
       'software':fm['software'],'model':model,'implementation':pins,'tasks':tasks,'tolerances':TOL,
       'source_full':record(full),'source_core':record(core),'numerical_reference':fm['numerical_reference'],
       'run_inventory':{'short_energy_gradient_calls':34,'new_DFT_calls':0},
       'density_coefficients':None,'charge_response_status':'not_part_of_this_component',
       'relaxation_correction_kcal_mol':None,'calibrated_class':None}
    for t in tasks:t['cache_key']=cache_key({'task':t,'model':model,'software':m['software'],'implementation':pins})
    write_new(out/'manifest.json',m)
    return validate(out/'manifest.json')


def validate(manifest):
    m=read_json(manifest)
    if m['schema_version']!=SCHEMA or m['tolerances']!=TOL:raise InvalidArtifact('short protocol/tolerances changed')
    sources={label:actual_collection(verify(m['source_'+label])) for label in ('full','core')}
    fc,fm=sources['full'];cc,cm=sources['core']
    expected=copy.deepcopy(fm['model']);expected.update(energy_component=COMPONENT,short_adapter=ADAPTER,
                                                     preparation_policy='exact_saved_full_and_GGR_curvature_v1')
    if m['model']!=expected or m['software']!=fm['software'] or numerical_parent_gate(m)['status']!='pass':
        raise InvalidArtifact('short model/software differs from validated parent')
    software=read_json(verify(m['software']))
    for ref in [m['agreement'],*m['implementation'].values(),software['python'],software['requirements'],
                *read_json(verify(software['backend_source_inventory']))['files']]:verify(ref)
    expected_ids={prefix+'__'+t['task_id'] for prefix,(_,source) in sources.items() for t in source['tasks']}
    if len(m['tasks'])!=34 or {t['task_id'] for t in m['tasks']}!=expected_ids:
        raise InvalidArtifact('short task inventory changed')
    for t in m['tasks']:
        prefix,source_id=t['task_id'].split('__',1);c,sm=sources[prefix]
        old=next(t for t in sm['tasks'] if t['task_id']==source_id)
        for key,value in old.items():
            if key not in ('task_id','cache_key') and t.get(key)!=value:
                raise InvalidArtifact('short physical input differs from saved full-forward input')
        if (t['energy_component']!=COMPONENT or t['source_collection']!=m['source_'+prefix] or
                t['reference_component_energy_eV']!=c['rows'][source_id]['energy_components_eV'][COMPONENT]):
            raise InvalidArtifact('short source/component reference changed')
        check_atoms(xyz(verify(t['xyz'])),t['charge'])
        base={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':base,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('scientific cache key mismatch')
    return {'status':'pass','tasks':34,'manifest':record(manifest)}


def collect_short(manifest):
    mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[];checks=[]
    for t in m['tasks']:
        valid=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            if r is not None:valid.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,
                             'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        r=valid[-1] if valid else {'status':'unavailable','energy_eV':None};rows[t['task_id']]=r
        error=None if r['status']!='computed' else r['energy_eV']-t['reference_component_energy_eV']
        checks.append({'name':t['task_id']+'_exact_component','energy_error_eV':error,
                       'pass':error is not None and abs(error)<=TOL['energy_eV']})
    complete=all(r['status']=='computed' for r in rows.values())
    if complete:
        for t in m['tasks']:
            if t.get('variant')!='rotate':continue
            primary=rows[t['task_id'].replace('_rotate','_primary')];r=rows[t['task_id']]
            energy_error=abs(r['energy_eV']-primary['energy_eV'])*EV_TO_KCAL
            a=np.load(verify(primary['forces']));b=np.load(verify(r['forces']))@np.array(t['rotation_matrix'])
            force_error=float(np.max(np.abs(a-b)))
            checks.append({'name':t['task_id']+'_rigid','energy_error_kcal_mol':energy_error,
                           'force_max_error_eV_A':force_error,'pass':energy_error<=TOL['energy_kcal_mol'] and force_error<=TOL['force_max_eV_A']})
        _,cm=actual_collection(verify(m['source_core']));_,dm=dft_source(verify(cm['DFT_collection']))
        for d in dm['directions']:
            rep,coord,h=d['representation'],d['coordinate'],d['amplitude'];jac=np.array(d['qm_jacobian'])
            values={}
            for metal in ('La','Ca'):
                minus,center,plus=[rows[f'core__{rep}_{name}_{metal}'] for name in (coord+'_minus','center',coord+'_plus')]
                odd=(plus['energy_eV']-minus['energy_eV'])*EV_TO_KCAL/2
                projected=-float(np.sum(np.load(verify(center['forces']))*jac))*EV_TO_KCAL
                values[metal]={'odd':odd,'predicted':h*projected}
            values['R']={k:values['Ca'][k]-values['La'][k] for k in ('odd','predicted')}
            for metal,v in values.items():
                err=v['odd']-v['predicted'];tol=max(TOL['odd_absolute_kcal_mol'],TOL['odd_relative']*abs(v['predicted']))
                checks.append({'name':f'{rep}_{coord}_{metal}_short_gradient','error_kcal_mol':err,
                               'tolerance_kcal_mol':tol,'pass':abs(err)<=tol})
    return {'status':'complete' if complete else 'incomplete','protocol_id':m['protocol_id'],'manifest':record(mp),
            'rows':rows,'attempts':attempts,'checks':checks,'numerical_checks_pass':complete and all(c['pass'] for c in checks),
            'energy_definition':'trained_local_interaction_energy_only','density_coefficients':None,
            'charge_response_status':'not_part_of_this_component','relaxation_correction_kcal_mol':None,'calibrated_class':None}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('full','core','agreement','output'):p.add_argument('--'+name,required=True)
    a=p.parse_args();print(json.dumps(prepare(**vars(a)),indent=2))
