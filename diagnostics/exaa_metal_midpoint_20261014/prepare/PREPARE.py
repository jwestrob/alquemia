#!/usr/bin/env python3
"""Metadata freeze and Slurm-only, single-fraction ion midpoint preparation."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import socket
import subprocess
import sys
import time
import traceback

PROTOCOL = 'ExaA_fixed_half_metal_displacement_saved_two_donor_geometries_v1'
INVENTORY_SHA = '757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343'
METALS = ('Ca','La')
CANDIDATES = ('adaptive_Ca','adaptive_La')
THREADS = ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','OPENMM_CPU_THREADS')


def pin(path):
    p = Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def check(p):
    if pin(p['path']) != p: raise ValueError('Changed pin: '+p['path'])
    return Path(p['path'])


def read(p):
    return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())


def write(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def freeze(root,output):
    root=Path(root).resolve();out=Path(output).resolve()
    inventory=pin(root/'workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json')
    if inventory['sha256'] != INVENTORY_SHA:raise ValueError('Wrong frozen source inventory')
    previous=pin(root/'workspaces/exaa_metal_position_20261014/preparation_v1/execution/PREPARATION.json')
    prior=read(previous);old=read(prior['manifest'])
    if prior['status']!='prepared':raise ValueError('Actual full-transplant preparation unavailable')
    sources=[]
    for sample in (0,3):
        cid=f'q9z4j7-pqq-la_model__conditioned_2Ca_v1__seed-1_sample-{sample}'
        src,=[r for r in read(inventory)['rows'] if r['case_id']==cid]
        prev,=[r for r in prior['rows'] if r['source_case_id']==cid]
        proposals={z:read(src['native_proposals'][z]) for z in METALS}
        if proposals['Ca']['manifest'] != proposals['La']['manifest']:raise ValueError('Original proposal manifests differ')
        sources.append({'sample':sample,'source_case_id':cid,'role':src['role'],
                        'record':src,'full_transplant_row':prev,
                        'proposal_manifest':proposals['Ca']['manifest']})
    out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir();pins={}
    for name in ('PREPARE.py','RUN.sbatch'):
        dst=impl/name;shutil.copyfile(Path(__file__).with_name(name),dst);pins[name]=pin(dst)
    m={'protocol':PROTOCOL,'inventory':inventory,'previous_full_transplant_preparation':previous,
       'previous_preparation_manifest':prior['manifest'],'sources':sources,
       'plan':pin(root/'diagnostics/exaa_metal_midpoint_20261014/PLAN.md'),
       'implementation':pins,'source_snapshot':old['source_snapshot'],'source_modules':old['source_modules'],
       'python':old['python'],'original_worker_settings':old['original_worker_settings'],
       'fraction':0.5,'candidate_order':list(CANDIDATES),'declared_sources':2,'declared_new_cells':8,
       'settings':{'kinematic_replay_A':1e-8,'maximum_metal_displacement_A':0.8,
                   'midpoint_formula':'original + fraction*(actual_full_transplant - original)',
                   'serialization':'.17g with exact float roundtrip'},
       'output':str(out/'execution'),'new_energy_calls':0,'new_optimizations':0,
       'new_alignment':False,'scientific_execution':'unrun','modeled_not_observed':True}
    write(out/'manifest.json',m);print(json.dumps(pin(out/'manifest.json')))


def placement():
    job=os.environ.get('SLURM_JOB_ID')
    if not job or not job.isdigit():raise RuntimeError('Preparation requires Slurm')
    text=subprocess.check_output(['scontrol','show','job','-o',job],text=True)
    d=dict(x.split('=',1) for x in text.split() if '=' in x)
    if (d.get('JobState'),d.get('NumNodes'),d.get('NumCPUs'))!=('RUNNING','1','1'):
        raise RuntimeError('Expected one CPU on one running node')
    if int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=1 or any(os.environ.get(k)!='1' for k in THREADS):
        raise RuntimeError('Expected one numerical thread')
    return {'job_id':job,'hostname':socket.gethostname(),'scontrol':text,'affinity':sorted(os.sched_getaffinity(0)),'allocated_CPUs':1}


def execute(manifest,sha):
    place=placement()  # Before numerical imports.
    mp={'path':str(Path(manifest).resolve()),'sha256':sha};m=read(mp)
    if m['protocol']!=PROTOCOL or m['fraction']!=0.5 or pin(__file__)!=m['implementation']['PREPARE.py']:
        raise ValueError('Wrong fixed midpoint protocol')
    out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);began=time.monotonic();cpu0=time.process_time()
    receipt={'status':'running','manifest':mp,'placement':place,'new_energy_calls':0,'new_optimizations':0}
    rows=[{'case_id':s['source_case_id'],'source_case_id':s['source_case_id'],'sample':s['sample'],'role':s['role'],
           'status':'unavailable','reason':'preparation not completed',
           'matrix':{z:{q:{'status':'unavailable','reason':'preparation not completed'} for q in CANDIDATES} for z in METALS}}
          for s in m['sources']]
    write(out/'START.json',receipt)
    try:
        check(m['inventory']);check(m['previous_full_transplant_preparation']);check(m['previous_preparation_manifest']);check(m['plan']);check(m['python'])
        if pin(sys.executable)!=m['python']:raise ValueError('Numerical interpreter differs')
        for p in m['implementation'].values():check(p)
        for p in m['source_modules'].values():check(p)
        sys.path.insert(0,m['source_snapshot'])
        import numpy as np
        import adaptive_metal_proposals as metal_geometry
        import mace_site_kinematics as kinematic_module
        import mace_hybrid as atomic_module
        import affordable_common as common
        for module in (metal_geometry,kinematic_module,atomic_module,common):
            if pin(module.__file__)!=m['source_modules'][Path(module.__file__).name]:raise ValueError('Wrong frozen helper import')
        if metal_geometry.SETTINGS['final_physical_tolerance_A']!=m['original_worker_settings']['final_physical_tolerance_A']:
            raise ValueError('Physical tolerance changed')
        inv=read(m['inventory']);prior=read(m['previous_full_transplant_preparation'])
        for s,row in zip(m['sources'],rows):
            src=s['record'];old=s['full_transplant_row'];cid=s['source_case_id'];directory=out/f"sample{s['sample']}";directory.mkdir()
            admission={'status':'unavailable','pass':False,'source_case_id':cid,'fraction':0.5,
                       'preserved_original_displacement_origin':True,'candidate_checks':{},'new_alignment':False}
            try:
                exact,=[x for x in inv['rows'] if x['case_id']==cid];prev,=[x for x in prior['rows'] if x['source_case_id']==cid]
                if src!=exact or old!=prev or old['status']!='prepared':raise ValueError('Saved source/transplant identity differs')
                for k in ('maps','state_signature','source_preparation'):
                    if old[k]!=src[k]:raise ValueError('Full-transplant source/state changed')
                for k in ('state_signature','source_preparation','prepared_admission','full_source','spectator_admission'):check(src[k])
                maps={z:read(src['maps'][z]) for z in METALS};state=read(src['state_signature']);a=read(src['prepared_admission'])
                if maps['Ca']!=maps['La'] or state!=src['state']:raise ValueError('Paired maps/state differs')
                if state['endpoint_charges']!={'Ca':-3,'La':-2} or state['multiplicities']!={'Ca':1,'La':1} or state['water_inventory']:
                    raise ValueError('Original dry chemical state differs')
                if a['status']!='prepared' or not a['mapping_eligible'] or not a['paired_maps_equal'] or (a['atom_count'],a['cap_count'])!=(192,14):
                    raise ValueError('Original preparation admission differs')
                if not read(old['geometry_admission'])['pass']:raise ValueError('Full transplant geometry was not admitted')
                kin=kinematic_module.Kinematics(maps['Ca']['context']);ids=[x['id'] for x in kin.modes]
                midx=[ids.index('metal_'+axis) for axis in 'xyz']
                for k,i in enumerate(midx):
                    mode=kin.modes[i]
                    if mode['kind']!='metal_translation' or mode['unit']!='angstrom' or mode['moving_indices']!=[kin.metal] or not np.array_equal(mode['axis'],np.eye(3)[k]):
                        raise ValueError('Unsupported metal translation modes')
                origin_atoms=common.xyz(check(src['matrix']['Ca']['origin']['xyz']))
                full_atoms=common.xyz(check(old['matrix']['Ca']['origin']['xyz']))
                if origin_atoms[1:]!=full_atoms[1:]:raise ValueError('Prior transplant moved nonmetal atoms')
                original=np.asarray(origin_atoms[0][1:]);full=np.asarray(full_atoms[0][1:]);midpoint=original+0.5*(full-original);delta=midpoint-original
                displacement=float(np.linalg.norm(delta))
                if displacement>0.8+metal_geometry.SETTINGS['final_physical_tolerance_A']:raise ValueError('Midpoint exceeds unchanged displacement bound')
                admission.update(original_metal_A=original.tolist(),actual_full_transplant_metal_A=full.tolist(),
                                 midpoint_metal_A=midpoint.tolist(),metal_displacement_A=displacement,
                                 original_origin_xyz=src['matrix']['Ca']['origin']['xyz'],
                                 actual_full_transplant_origin_xyz=old['matrix']['Ca']['origin']['xyz'])
                row.update({k:src[k] for k in ('selection_id','maps','state_signature','source_preparation','prepared_admission','full_source','spectator_admission')})
                row.update(fraction=0.5,original_metal_A=original.tolist(),actual_full_transplant_metal_A=full.tolist(),
                           midpoint_metal_A=midpoint.tolist(),modeled_not_observed=True,
                           original_origin_xyz=src['matrix']['Ca']['origin']['xyz'],
                           actual_full_transplant_origin_xyz=old['matrix']['Ca']['origin']['xyz'],
                           prior_full_transplant_admission=old['geometry_admission'])
                pm=read(s['proposal_manifest']);tasks={t['metal']:t for t in pm['tasks']}
                for name in CANDIDATES:
                    owner=name.rsplit('_',1)[1];proposal=read(src['native_proposals'][owner]);task=copy.deepcopy(tasks[owner])
                    if proposal['status']!='proposal_available' or proposal['manifest']!=s['proposal_manifest'] or not proposal['proposal']['geometry_checks']['pass']:
                        raise ValueError('Original donor geometry unavailable')
                    q=np.asarray(proposal['proposal']['full_q']);task['active_indices']=midx+task['active_indices']
                    if len(task['active_indices'])!=7 or len(set(task['active_indices']))!=7 or not np.array_equal(q[midx],np.zeros(3)):
                        raise ValueError('Original active donor subspace differs')
                    if not np.array_equal(metal_geometry.full_q(task,q[task['active_indices']]),q):raise ValueError('Saved q outside original donor subspace')
                    oldca=common.xyz(check(src['matrix']['Ca'][name]['xyz']));oldla=common.xyz(check(src['matrix']['La'][name]['xyz']))
                    if oldca[1:]!=oldla[1:] or oldca[0][1:]!=oldla[0][1:] or not np.array_equal(oldca[0][1:],original):
                        raise ValueError('Original paired coordinates/metal differ')
                    old_error=float(np.max(np.abs(kin.evaluate(q)[1]-np.asarray([x[1:] for x in oldca]))))
                    if old_error>m['settings']['kinematic_replay_A']:raise ValueError('Original donor geometry does not replay')
                    oldcheck=metal_geometry.final_geometry(kin,task,q[task['active_indices']],[x[0] for x in oldca])
                    q=q.copy();q[midx]=delta
                    physical=metal_geometry.final_geometry(kin,task,q[task['active_indices']],[x[0] for x in oldca])
                    expected=np.asarray([x[1:] for x in oldca]);expected[0]=midpoint
                    error=float(np.max(np.abs(kin.evaluate(q)[1]-expected)))
                    if error>m['settings']['kinematic_replay_A']:raise ValueError('Midpoint replay changes nonmetal geometry')
                    gate={'status':'prepared','pass':True,'original_donor_candidate':name,
                          'original_donor_proposal':src['native_proposals'][owner],'original_geometry':oldcheck,
                          'physical_geometry':physical,'full_q':q.tolist(),'original_replay_max_A':old_error,
                          'rendered_vs_exact_saved_max_A':error,'fraction':0.5}
                    admission['candidate_checks'][name]=gate
                    for z in METALS:
                        cell=src['matrix'][z][name];xp=directory/f'{z}_{name}.xyz';original_path=check(cell['xyz'])
                        if (cell['atom_count'],cell['multiplicity'],cell['GFN2_valence_electrons'])!=(192,1,588) or cell['charge']!=state['endpoint_charges'][z]:
                            raise ValueError('Original endpoint state differs')
                        check(cell['MACE']);check(cell['MACE_request']);check(cell['native_frozen_ALPB_result'])
                        lines=original_path.read_text().splitlines(keepends=True);original_lines=list(lines)
                        lines[2]=z+' '+' '.join(format(float(v),'.17g') for v in midpoint)+'\n'
                        with xp.open('x') as f:f.write(''.join(lines))
                        actual=common.xyz(xp);atomic=atomic_module.check_atoms(actual,cell['charge'])
                        if xp.read_text().splitlines(keepends=True)[3:]!=original_lines[3:] or actual[1:]!=(oldca if z=='Ca' else oldla)[1:]:
                            raise ValueError('Nonmetal lines/coordinates changed')
                        if not np.array_equal(actual[0][1:],midpoint):raise ValueError('Midpoint serialization changed')
                        row['matrix'][z][name]={'status':'prepared','xyz':pin(xp),'charge':cell['charge'],'multiplicity':1,
                            'atom_count':192,'valence_electrons':588,'atomic_check':atomic,'original_donor_candidate':name,
                            'original_xyz':cell['xyz'],'original_MACE':cell['MACE'],'original_MACE_eV':cell['MACE_eV'],
                            'original_frozen_ALPB':cell['native_frozen_ALPB_result'],'geometry_checks':gate}
                    ca=common.xyz(check(row['matrix']['Ca'][name]['xyz']));la=common.xyz(check(row['matrix']['La'][name]['xyz']))
                    if ca[1:]!=la[1:] or ca[0][1:]!=la[0][1:]:raise ValueError('New paired coordinates differ')
                admission.update(status='prepared',**{'pass':True},nonmetal_lines_unchanged=True,paired_geometry_identical=True)
                row.update(status='prepared',reason=None)
            except Exception as exc:
                admission.update(reason=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc());row['reason']=admission['reason']
                for z in METALS:
                    for name in CANDIDATES:row['matrix'][z][name].update(status='unavailable',reason=admission['reason'])
            write(directory/'GEOMETRY_ADMISSION.json',admission);row['geometry_admission']=pin(directory/'GEOMETRY_ADMISSION.json')
            for z in METALS:
                for name in CANDIDATES:row['matrix'][z][name]['geometry_admission']=row['geometry_admission']
        receipt['status']='complete'
    except BaseException as exc:
        receipt.update(status='failed',error=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc());raise
    finally:
        receipt.update(wall_seconds=time.monotonic()-began,process_CPU_seconds=time.process_time()-cpu0,
                       peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(out/'EXECUTION.json',receipt)
        result={'protocol':PROTOCOL,'manifest':mp,'inventory':m['inventory'],
                'previous_full_transplant_preparation':m['previous_full_transplant_preparation'],
                'terminal':True,'status':'prepared' if all(r['status']=='prepared' for r in rows) else 'unavailable',
                'rows':rows,'declared_sources':2,'declared_new_cells':8,'prepared_sources':sum(r['status']=='prepared' for r in rows),
                'execution':pin(out/'EXECUTION.json'),'fraction':0.5,'new_energy_calls':0,'new_optimizations':0,
                'new_alignment':False,'modeled_not_observed':True,'production_changed':False}
        write(out/'PREPARATION.json',result)
    print(json.dumps({'preparation':pin(out/'PREPARATION.json'),'prepared_sources':result['prepared_sources']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('operation',choices=('freeze','execute'))
    p.add_argument('--root');p.add_argument('--output');p.add_argument('--manifest');p.add_argument('--manifest-sha')
    a=p.parse_args()
    if a.operation=='freeze':freeze(a.root,a.output)
    else:execute(a.manifest,a.manifest_sha)
