"""Bounded real-coordinate qualification of analytic masked-score derivatives."""
from __future__ import annotations
import argparse
import copy
import csv
import json
import math
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL, accepted_attempt, check_atoms, rotation, write_xyz
from mace_omol_ablation import ADAPTER as MASK, COMPONENT
from mace_omol_ablation_run import PROTOCOL, SEMANTICS, descriptor_model
from mace_omol_gradients import ADAPTER

STAGES = ('masked_gradient_core','masked_gradient_full')
CONFIG = {'id':'masked_descriptor_gradient_check_v1','step_A':.01,
          'coordinate_replay_tolerance_A':1e-12,
          'energy_tolerance_model_kcal':.01,'gradient_tolerance_eV_A':.001,
          'odd_absolute_tolerance_model_kcal':.01,'odd_relative_tolerance':.01,
          'checkpoint_adapter':ADAPTER,'chunk_size':1024,'checkpoint_early_stop':False}
REPORT_REPLAY_TOLERANCE=1e-10  # Host arithmetic only; far below every scientific gate.


def geometry(task):
    source = xyz(verify(task['source_xyz'])); points=np.array([a[1:] for a in source])
    if task['gradient_variant']=='rotate': points=points@rotation().T
    elif task['gradient_variant'] in ('positive','negative'):
        sign=1 if task['gradient_variant']=='positive' else -1
        points[task['metal_index']]+=sign*CONFIG['step_A']*np.array(task['direction'])
    return [(a[0],*p) for a,p in zip(source,points)]


@cached_file_checks
def expected(development_collection, stage):
    from mace_omol_prepared import qualification
    saved,parent=qualification(development_collection)
    if stage not in STAGES: raise InvalidArtifact('unsupported derivative stage')
    tasks=[]
    for metal in ('Ca','La'):
        source_id=(f'bridge_1H4I_{metal}_mask_native' if stage==STAGES[0]
                   else f'GGR_1GLG_{metal}_bound_primary_mask_batched')
        source=next(t for t in parent['tasks'] if t['task_id']==source_id)
        coords=xyz(verify(source['xyz'])); index=source['metal_index']; pos=np.array([a[1:] for a in coords])
        if stage==STAGES[0]:
            oxygen=min((i for i,a in enumerate(coords) if a[0]=='O'),key=lambda i:(np.linalg.norm(pos[i]-pos[index]),i))
            selector={'kind':'nearest_source_oxygen','source_index':oxygen}
        else:
            p=read_json(verify(source['preparation']))
            oxygen=next(i for i,a in enumerate(p['physical_atoms']) if a['id']=='A/140//O')
            selector={'kind':'source_atom_id','id':'A/140//O','source_index':oxygen}
        direction=pos[oxygen]-pos[index]; direction/=np.linalg.norm(direction)
        variants=[('checkpointed','center'),('checkpointed','positive'),('checkpointed','negative')]
        if stage==STAGES[0]:variants=[('native','center'),('checkpointed','center'),('checkpointed','rotate'),*variants[1:]]
        for backend, variant in variants:
            task={k:copy.deepcopy(source[k]) for k in ('case_id','metal','kind','metal_index','charge',
                  'spin_multiplicity','state','preparation','assembly','microstate','explicit_waters','evidence') if k in source}
            task.update(task_id=f'{metal}_{backend}_{variant}',source_xyz=source['xyz'],
                        source_task={'manifest':saved['manifest'],'task_id':source_id},
                        descriptor_gradient_experiment=CONFIG['id'],derivative_backend=backend,
                        gradient_variant=variant,direction=direction.tolist(),direction_selector=selector,
                        energy_only=variant in ('positive','negative'),energy_component=COMPONENT,
                        charge_feature_adapter=MASK,output_semantics=SEMANTICS,capture_native_readout=True,
                        archived_center_energy_eV=saved['rows'][source_id]['energy_eV'])
            tasks.append(task)
    return tasks,parent


def model(software):
    m=descriptor_model(software);m['gradient_execution_adapter']=CONFIG
    return m


@cached_file_checks
def prepare(development_collection, agreement, output, stage, core_report=None):
    from mace_omol import common,seal
    if stage==STAGES[1] and core_report is None:
        raise InvalidArtifact('passing actual core-gradient qualification required')
    tasks,parent=expected(development_collection,stage)
    _,out,m=common(verify(parent['inventory']),verify(parent['software']),agreement,output,stage)
    m.update(protocol_id=PROTOCOL,model=model(verify(parent['software'])),gradient_settings=CONFIG,
             development_collection=record(development_collection),tasks=tasks,
             core_report=record(core_report) if core_report else None,output_semantics=SEMANTICS)
    if stage==STAGES[1]: check_parent(m)
    elif core_report is not None: raise InvalidArtifact('core stage cannot borrow a qualification')
    for t in tasks:
        path=out/(t['task_id']+'.xyz'); write_xyz(path,geometry(t));t['xyz']=record(path)
    return seal(out,m)


def check_parent(manifest):
    p=read_json(verify(manifest['core_report']));mp=verify(p['manifest']);m=read_json(mp)
    if m['stage']!=STAGES[0] or not p['numerical_gate_pass'] or p['status']!='complete':
        raise InvalidArtifact('passing actual core-gradient qualification required')
    actual=collect(mp)
    # Scientific artifacts, settings, success decisions and receipts stay exact.
    # BLAS reductions of those same arrays need not be bit-identical on two CPUs.
    for key in actual:
        if key not in ('checks','projections') and p[key]!=actual[key]:
            raise InvalidArtifact('core report differs from actual receipts: '+key)
    differences=[]
    def derived(a,b,path):
        if isinstance(a,dict) and isinstance(b,dict) and a.keys()==b.keys():
            for key in a:derived(a[key],b[key],path+'/'+key)
        elif isinstance(a,list) and isinstance(b,list) and len(a)==len(b):
            for i,(x,y) in enumerate(zip(a,b)):derived(x,y,path+'/'+str(i))
        elif isinstance(a,float) and isinstance(b,float):
            error=abs(a-b)
            if not math.isfinite(error) or error>REPORT_REPLAY_TOLERANCE:
                raise InvalidArtifact(f'core derived report differs at {path}: {a} versus {b}')
            if error:differences.append({'field':path,'absolute_difference':error})
        elif a!=b:raise InvalidArtifact('core derived report type/value differs at '+path)
    if len(p['checks'])!=len(actual['checks']):raise InvalidArtifact('core check inventory differs')
    for i,(old,new) in enumerate(zip(p['checks'],actual['checks'])):
        if {k:v for k,v in old.items() if k!='error'}!={k:v for k,v in new.items() if k!='error'}:
            raise InvalidArtifact('core criterion or pass/fail decision changed')
        derived(old['error'],new['error'],f'checks/{i}/error')
    derived(p['projections'],actual['projections'],'projections')
    if differences:
        print(json.dumps({'event':'derived_report_replay_roundoff','tolerance':REPORT_REPLAY_TOLERANCE,
                          'differences':differences}),flush=True)
    for name in ('mace_omol_gradients.py','mace_omol_gradient_worker.py'):
        if manifest['implementation'][name]['sha256']!=m['implementation'][name]['sha256']:
            raise InvalidArtifact('gradient execution changed after qualification')


@cached_file_checks
def validate(manifest):
    m=read_json(manifest);tasks,parent=expected(verify(m['development_collection']),m['stage'])
    if (m['protocol_id']!=PROTOCOL or m['gradient_settings']!=CONFIG or m['reused']
            or m['model']!=model(verify(m['software'])) or m['software']!=parent['software']
            or m['inventory']!=parent['inventory'] or len(m['tasks'])!=len(tasks)):
        raise InvalidArtifact('derivative experiment inputs or model changed')
    for pin in [m['agreement'],*m['implementation'].values()]: verify(pin)
    if m['stage']==STAGES[1]: check_parent(m)
    elif m['core_report'] is not None: raise InvalidArtifact('unexpected core qualification reference')
    maximum_coordinate_error=0.;coordinate_roundoff=[]
    for actual,wanted in zip(m['tasks'],tasks):
        if {k:v for k,v in actual.items() if k not in ('xyz','cache_key')}!=wanted:
            raise InvalidArtifact('gradient task differs from fixed inventory')
        coords=xyz(verify(actual['xyz']))
        replay=geometry(wanted)
        if len(coords)!=len(replay) or [a[0] for a in coords]!=[a[0] for a in replay]:
            raise InvalidArtifact('gradient atom inventory differs')
        error=float(np.max(np.abs(np.array([a[1:] for a in coords])-np.array([a[1:] for a in replay]))))
        maximum_coordinate_error=max(maximum_coordinate_error,error)
        if error:coordinate_roundoff.append({'task_id':actual['task_id'],'max_error_A':error})
        if error>CONFIG['coordinate_replay_tolerance_A'] or check_atoms(coords,actual['charge'])!=actual['state']:
            raise InvalidArtifact(f"gradient coordinates/charge differ for {actual['task_id']}: max_error_A={error}")
        payload={k:v for k,v in actual.items() if k!='cache_key'}
        if actual['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('derivative cache identity differs')
    if coordinate_roundoff:
        print(json.dumps({'event':'coordinate_replay_roundoff','tolerance_A':CONFIG['coordinate_replay_tolerance_A'],
                          'rows':coordinate_roundoff}),flush=True)
    return {'status':'pass','tasks':len(tasks),'manifest':record(manifest),
            'maximum_coordinate_replay_error_A':maximum_coordinate_error}


@cached_file_checks
def collect(manifest):
    mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[];checks=[];projections={}
    for task in m['tasks']:
        found=[]
        for path in sorted((mp.parent/'execution'/task['task_id']).glob('attempt_*')):
            row=accepted_attempt(path,task,mp)
            if row is not None: found.append(row)
            attempts.append({'task_id':task['task_id'],'path':str(path),'accepted':row is not None,
                             'receipt':record(path/'receipt.json') if (path/'receipt.json').exists() else None})
        rows[task['task_id']]=found[-1] if found else {'status':'unavailable','energy_eV':None,'forces':None}
    complete=all(r['status']=='computed' for r in rows.values())
    def check(name,error,tolerance):
        error=float(error);tolerance=float(tolerance)
        checks.append({'name':name,'error':error,'tolerance':tolerance,'pass':abs(error)<=tolerance})
    if complete:
        for name,row in rows.items():check(name+'_readout_model_kcal',row['native_readout']['component_sum_error_kcal_mol'],.01)
        for metal in ('Ca','La'):
            center=rows[metal+'_checkpointed_center'];forces=np.load(verify(center['forces']))
            task=next(t for t in m['tasks'] if t['task_id']==metal+'_checkpointed_center')
            check(metal+'_archived_center_model_kcal',(center['energy_eV']-task['archived_center_energy_eV'])*EV_TO_KCAL,.01)
            if m['stage']==STAGES[0]:
                native=rows[metal+'_native_center'];rotated=rows[metal+'_checkpointed_rotate']
                check(metal+'_native_energy_model_kcal',(center['energy_eV']-native['energy_eV'])*EV_TO_KCAL,.01)
                check(metal+'_native_gradient_eV_A',np.max(np.abs(forces-np.load(verify(native['forces'])))),.001)
                check(metal+'_rotated_energy_model_kcal',(rotated['energy_eV']-center['energy_eV'])*EV_TO_KCAL,.01)
                check(metal+'_rotated_gradient_eV_A',np.max(np.abs(np.load(verify(rotated['forces']))@rotation()-forces)),.001)
            gradient=-forces;g=float(gradient[task['metal_index']]@np.array(task['direction']))
            odd=.5*(rows[metal+'_checkpointed_positive']['energy_eV']-rows[metal+'_checkpointed_negative']['energy_eV'])*EV_TO_KCAL
            prediction=CONFIG['step_A']*g*EV_TO_KCAL
            check(metal+'_odd_model_kcal',odd-prediction,max(.01,.01*abs(prediction)))
            projections[metal]={'gradient_eV_equivalent_per_A':g,'odd_model_kcal':odd,'predicted_odd_model_kcal':prediction,
                                'max_atom_gradient_norm_eV_A':float(np.linalg.norm(gradient,axis=1).max())}
        odd=projections['Ca']['odd_model_kcal']-projections['La']['odd_model_kcal']
        prediction=projections['Ca']['predicted_odd_model_kcal']-projections['La']['predicted_odd_model_kcal']
        check('R_odd_model_kcal',odd-prediction,max(.01,.01*abs(prediction)))
        projections['R']={'gradient_model_kcal_per_A':(projections['Ca']['gradient_eV_equivalent_per_A']-projections['La']['gradient_eV_equivalent_per_A'])*EV_TO_KCAL,
                          'odd_model_kcal':odd,'predicted_odd_model_kcal':prediction}
    return {'status':'complete' if complete else 'incomplete','manifest':record(mp),'protocol_id':PROTOCOL,
            'derivative_experiment':CONFIG,'rows':rows,'attempts':attempts,'checks':checks,
            'numerical_gate_pass':complete and all(c['pass'] for c in checks),'projections':projections,
            'physical_force_validation_claimed':False,'relaxation_correction_kcal_mol':None,
            'baseline_changed':False,'biological_score_changes_applied':False}


@cached_file_checks
def report(manifest,output):
    validate(manifest);r=collect(manifest);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    if r['status']=='complete':
        ca,la=(np.load(verify(r['rows'][m+'_checkpointed_center']['forces'])) for m in ('Ca','La'))
        path=out/'paired_descriptor_gradient_model_kcal_per_A.npy';np.save(path,(-ca+la)*EV_TO_KCAL)
        r['paired_gradient']=record(path)
        task=next(t for t in read_json(manifest)['tasks'] if t['task_id']=='Ca_checkpointed_center')
        r['gradient_coordinate_source']=task['source_xyz'];r['preparation']=task.get('preparation')
        if r['preparation'] is not None:
            atoms=read_json(verify(r['preparation']))['physical_atoms']
            coords=xyz(verify(task['source_xyz']))
            if (len(atoms)!=len(coords) or len(atoms)!=len(ca)
                    or len({a['id'] for a in atoms})!=len(atoms)):
                raise InvalidArtifact('gradient source mapping is not bijective')
            for i,(atom,coordinate) in enumerate(zip(atoms,coords)):
                element=task['metal'] if atom['kind']=='selected_metal' else atom['element']
                if (element!=coordinate[0] or np.max(np.abs(np.array(atom['xyz_A'])-coordinate[1:]))>1e-12):
                    raise InvalidArtifact('gradient source coordinate differs at '+str(i))
            path=out/'source_mapped_descriptor_gradients.tsv'
            with path.open('x',newline='') as handle:
                writer=csv.writer(handle,delimiter='\t')
                writer.writerow(['index','source_id','kind','resname','element_Ca_endpoint',
                    'x_A','y_A','z_A',*[f'{metal}_d{axis}_model_kcal_per_A' for metal in ('Ca','La','R') for axis in ('x','y','z')]])
                for i,(atom,coordinate) in enumerate(zip(atoms,coords)):
                    writer.writerow([i,atom['id'],atom['kind'],atom.get('resname',''),coordinate[0],
                        *coordinate[1:],*(-ca[i]*EV_TO_KCAL),*(-la[i]*EV_TO_KCAL),*((-ca[i]+la[i])*EV_TO_KCAL)])
            r['source_mapped_gradients']=record(path)
    r['report_implementation']=record(__file__);write_new(out/'result.json',r)
    return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for key in ('development-collection','agreement','output'):a.add_argument('--'+key,required=True)
    a.add_argument('--stage',choices=STAGES,required=True);a.add_argument('--core-report')
    a=sub.add_parser('report');a.add_argument('--manifest',required=True);a.add_argument('--output',required=True)
    args=vars(p.parse_args());command=args.pop('command');r={'prepare':prepare,'report':report}[command](**args)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','attempts','checks')},indent=2))
