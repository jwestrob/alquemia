"""Native energy-only OMOL qualification and intact-chain coordination descriptor."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, rotation, accepted_attempt, check_atoms
from mace_omol_coordination import bound_source
from mace_global_benchmark import physical_preparations
from mace_global_prepare import CASES

PROTOCOL='mace_omol_intact_chain_matched_coordination_v1'
PREPARATION_SHA='63d6a44fbef3a1dc888696bba5d122ac3531fe12ceee4ee03a85194f768b0552'
EVALUATION={'forward':'native_ScaleShiftMACE','precision':'float64','compute_force':False,
            'torch_gradient_mode':'no_grad','native_atomic_reference':'single_omol_head',
            'forces_status':'not_requested','detached_offset_A':30.,'farther_offset_A':10.}
STAGES=('intact_core','intact_qualification','intact_benchmark')


def preparations(path):
    if record(path)['sha256']!=PREPARATION_SHA:
        raise InvalidArtifact('intact-chain physical source differs from frozen preparation')
    p,_=physical_preparations(path);result={}
    for r in p['rows']:
        c=read_json(verify(r['preparation']));indices=[i for i,a in enumerate(c['physical_atoms']) if a['element']=='M']
        if len(indices)!=1 or any('cap' in a['kind'] for a in c['physical_atoms']):
            raise InvalidArtifact('single selected metal and no synthetic carve caps required')
        result[r['case_id']]={'source':r['preparation'],'data':c,'metal_index':indices[0]}
    return result


def energy_task(task):
    t=copy.deepcopy(task);t.pop('cache_key',None)
    t.update(energy_only=True,capture_native_readout=True,requested_properties=['energy'])
    return t


def expected_tasks(stage, bound_tasks, physical):
    if stage=='intact_core':
        tasks=[]
        for t in bound_tasks:
            if t['case_id'] in ('1H4I','4MAE'):
                new=energy_task(t);new.update(task_id='bridge_'+t['task_id'],source_task_id=t['task_id'],
                                             source_xyz=t['xyz'],metal_index=0,position='bound')
                tasks.append(new)
        return tasks
    from mace_omol import COMPONENT
    tasks=[]
    for name in CASES:
        if stage=='intact_qualification' and name!='ALPHA_1F6S':continue
        prep=physical[name];c=prep['data']
        for position in ('bound','detached'):
            variants=('primary','repeat','rotate')+(('farther',) if position=='detached' else ()) if stage=='intact_qualification' else ('primary',)
            for metal in ('La','Ca'):
                e=c['endpoints'][metal]
                for variant in variants:
                    t={'task_id':f'{name}_{metal}_{position}_{variant}','case_id':name,'metal':metal,'kind':'full',
                       'position':position,'variant':variant,'metal_index':prep['metal_index'],'source_xyz':e['xyz'],
                       'preparation':prep['source'],'assembly':c['assembly'],'microstate':c['microstate'],
                       'explicit_waters':c['explicit_waters'],'evidence':c['evidence'],'energy_component':COMPONENT,
                       'require_isolated_metal':position=='detached',**{k:v for k,v in e.items() if k!='xyz'}}
                    tasks.append(energy_task(t))
    return tasks


def geometry(task):
    rows=xyz(verify(task['source_xyz']));i=task['metal_index']
    if rows[i][0]!=task['metal']:raise InvalidArtifact('physical selected-metal mapping differs')
    pos=np.array([r[1:] for r in rows]);others=np.delete(pos,i,axis=0)
    if task['position']=='detached':
        pos[i,0]=max(others[:,0])+EVALUATION['detached_offset_A']
        if task['variant']=='farther':pos[i,0]+=EVALUATION['farther_offset_A']
    if task['variant']=='rotate':pos=(pos-pos[i])@rotation().T+pos[i]
    if task['position']=='detached' and np.min(np.linalg.norm(np.delete(pos,i,axis=0)-pos[i],axis=1))<=6.:
        raise InvalidArtifact('detached metal still inside native edge cutoff')
    return [(row[0],*p) for row,p in zip(rows,pos)]


def qualified(path,stage):
    c=read_json(path);mp=verify(c['manifest']);m=read_json(mp)
    if m['stage']!=stage:raise InvalidArtifact('wrong intact qualification stage')
    validate(mp);actual=collect(mp)
    if c!=actual or not c['numerical_gate_pass']:
        raise InvalidArtifact('actual passing intact numerical qualification required')
    return c,m


def prepare(collection,preparation,agreement,output,core_qualification=None,full_qualification=None):
    from mace_omol import common,seal
    if full_qualification and not core_qualification:raise InvalidArtifact('full qualification requires its core bridge')
    _,parent,bt,_=bound_source(collection);physical=preparations(preparation)
    stage='intact_benchmark' if full_qualification else 'intact_qualification' if core_qualification else 'intact_core'
    _,out,m=common(verify(parent['inventory']),verify(parent['software']),agreement,output,stage)
    m.update(protocol_id=PROTOCOL,energy_evaluation=EVALUATION,source_collection=record(collection),
             physical_preparation=record(preparation))
    if core_qualification:
        qc,qm=qualified(core_qualification,'intact_core');m['core_qualification']=record(core_qualification)
        if (qm['source_collection']!=m['source_collection'] or qm['physical_preparation']!=m['physical_preparation']
                or qm['model']!=m['model'] or qm['software']!=m['software']):
            raise InvalidArtifact('core bridge belongs to another method/source')
    if full_qualification:
        qc,qm=qualified(full_qualification,'intact_qualification');m['full_qualification']=record(full_qualification)
        if qm['core_qualification']!=m['core_qualification']:
            raise InvalidArtifact('full qualification used another core bridge')
        m['reused']={t['task_id']:{'source_collection':record(full_qualification),'source_task_id':t['task_id']}
                     for t in qm['tasks'] if t['variant']=='primary'}
    wanted=[t for t in expected_tasks(stage,bt,physical) if t['task_id'] not in m['reused']]
    for t in wanted:
        if t['position']=='bound' and t['variant'] in ('primary','repeat'):
            t['xyz']=t['source_xyz']
        else:
            p=out/(t['task_id']+'.xyz');rows=geometry(t)
            with p.open('x') as f:
                f.write(f'{len(rows)}\nDerived intact-chain diagnostic; transformation pinned in manifest\n')
                for row in rows:f.write(row[0]+' '+' '.join(format(v,'.17g') for v in row[1:])+'\n')
            t['xyz']=record(p)
    m['tasks']=wanted
    return seal(out,m)


def validate(manifest):
    from mace_omol import SCHEMA,TOL,model
    m=read_json(manifest)
    if (m['schema_version']!=SCHEMA or m['protocol_id']!=PROTOCOL or m['stage'] not in STAGES
            or m['energy_evaluation']!=EVALUATION or m['tolerances']!=TOL):
        raise InvalidArtifact('intact energy-only model/acceptance changed')
    _,parent,bt,_=bound_source(verify(m['source_collection']));physical=preparations(verify(m['physical_preparation']))
    if (m['model']!=parent['model'] or m['software']!=parent['software'] or m['inventory']!=parent['inventory']
            or m['model']!=model(verify(m['software']))):
        raise InvalidArtifact('intact Hamiltonian differs from qualified OMOL')
    s=read_json(verify(m['software']))
    for pin in [m['agreement'],*m['implementation'].values(),s['python'],s['requirements'],s['checkpoint'],s['inspection'],
                s['download_receipt'],*read_json(verify(s['backend_source_inventory']))['files']]:verify(pin)
    if m['stage']!='intact_core':
        qc,qm=qualified(verify(m['core_qualification']),'intact_core')
        if (qm['source_collection']!=m['source_collection'] or qm['physical_preparation']!=m['physical_preparation']
                or qm['model']!=m['model'] or qm['software']!=m['software']):
            raise InvalidArtifact('core bridge source or method differs')
    if m['stage']=='intact_benchmark':
        qc,qm=qualified(verify(m['full_qualification']),'intact_qualification')
        if qm['core_qualification']!=m['core_qualification']:raise InvalidArtifact('full bridge source differs')
        reuse={t['task_id']:{'source_collection':m['full_qualification'],'source_task_id':t['task_id']}
               for t in qm['tasks'] if t['variant']=='primary'}
        if m['reused']!=reuse or len(reuse)!=4:raise InvalidArtifact('qualified intact reuse differs')
    elif m['reused']:raise InvalidArtifact('unexpected intact cache reuse')
    expected=[t for t in expected_tasks(m['stage'],bt,physical) if t['task_id'] not in m['reused']]
    if len(expected)!={'intact_core':4,'intact_qualification':14,'intact_benchmark':16}[m['stage']] or len(expected)!=len(m['tasks']):
        raise InvalidArtifact('finite intact task count changed')
    for t,e in zip(m['tasks'],expected):
        if any(t.get(k)!=v for k,v in e.items() if k!='xyz'):raise InvalidArtifact('intact task differs from declared source')
        actual=xyz(verify(t['xyz']));wanted=geometry(t)
        if [r[0] for r in actual]!=[r[0] for r in wanted] or np.max(np.abs(np.array([r[1:] for r in actual])-np.array([r[1:] for r in wanted])))>1e-12:
            raise InvalidArtifact('intact coordinates differ from physical inventory/declared transform')
        if t['position']=='bound' and t['variant'] in ('primary','repeat') and t['xyz']!=t['source_xyz']:
            raise InvalidArtifact('bound source bytes changed')
        if check_atoms(actual,t['charge'])!=t['state']:raise InvalidArtifact('state electron parity/count differs')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('intact scientific cache identity changed')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest)}


def collect(manifest):
    from mace_omol import TOL,UNAVAILABLE
    mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[];checks=[];reused={}
    for t in m['tasks']:
        found=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            if r is not None:found.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,
                             'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        rows[t['task_id']]=found[-1] if found else {'status':'unavailable','energy_eV':None,'forces':None}
    complete=all(r['status']=='computed' for r in rows.values())
    if complete:
        for key,r in rows.items():
            error=r['native_readout']['component_sum_error_kcal_mol']
            checks.append({'name':key+'_native_energy_accounting','error_kcal_mol':error,
                           'pass':abs(error)<=TOL['energy_kcal_mol'] and r['gradient_computation_enabled'] is False})
    if complete and m['stage']=='intact_core':
        _,_,_,bound=bound_source(verify(m['source_collection']));errors={}
        for t in m['tasks']:
            error=(rows[t['task_id']]['energy_eV']-bound[t['source_task_id']]['energy_eV'])*EV_TO_KCAL
            errors[t['source_task_id']]=error
            checks.append({'name':t['task_id']+'_force_producing_reference','error_kcal_mol':error,
                           'pass':abs(error)<=TOL['energy_kcal_mol']})
        for case in ('1H4I','4MAE'):
            error=errors[case+'_Ca']-errors[case+'_La']
            checks.append({'name':case+'_paired_R_bridge','error_kcal_mol':error,'pass':abs(error)<=TOL['energy_kcal_mol']})
    if complete and m['stage']=='intact_qualification':
        def value(metal,position,variant):return rows[f'ALPHA_1F6S_{metal}_{position}_{variant}']['energy_eV']
        errors={}
        for variant in ('repeat','rotate','farther'):
            for position in (('detached',) if variant=='farther' else ('bound','detached')):
                for metal in ('La','Ca'):
                    error=(value(metal,position,variant)-value(metal,position,'primary'))*EV_TO_KCAL
                    errors[(metal,position,variant)]=error
                    checks.append({'name':metal+'_'+position+'_'+variant,'error_kcal_mol':error,'pass':abs(error)<=TOL['energy_kcal_mol']})
                error=errors[('Ca',position,variant)]-errors[('La',position,variant)]
                checks.append({'name':'R_'+position+'_'+variant,'error_kcal_mol':error,'pass':abs(error)<=TOL['energy_kcal_mol']})
            if variant!='farther':
                error=(errors[('Ca','bound',variant)]-errors[('La','bound',variant)]
                       -errors[('Ca','detached',variant)]+errors[('La','detached',variant)])
                checks.append({'name':'R_coord_'+variant,'error_kcal_mol':error,'pass':abs(error)<=TOL['energy_kcal_mol']})
    if m['stage']=='intact_benchmark':
        qc,_=qualified(verify(m['full_qualification']),'intact_qualification')
        for key,ref in m['reused'].items():reused[key]={**ref,'result':qc['rows'][ref['source_task_id']]}
    return {'status':'complete' if complete else 'incomplete','manifest':record(mp),'rows':rows,'attempts':attempts,
            'reused_rows':reused,'checks':checks,'numerical_gate_pass':complete and bool(checks) and all(c['pass'] for c in checks),
            'force_validation_status':'not_requested_energy_only',**UNAVAILABLE}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('collection','preparation','agreement','output'):p.add_argument('--'+key,required=True)
    p.add_argument('--core-qualification');p.add_argument('--full-qualification')
    print(json.dumps(prepare(**vars(p.parse_args())),indent=2))
