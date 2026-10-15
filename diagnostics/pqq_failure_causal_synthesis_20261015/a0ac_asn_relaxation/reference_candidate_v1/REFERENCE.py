#!/usr/bin/env python3
"""Four fixed Ca4 selected energies, using the unchanged pinned native executor."""
import argparse
from contextlib import ExitStack
import datetime as dt
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import sys
import traceback

PROTOCOL = 'A0AC_Ca4_selected_old_new_native_r2SCAN3c_vacuum_v1'
SOURCE = 'a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4'
SELECTED = {'Ca': ['adaptive_Ca','relaxed_flip_Ca'], 'La': ['adaptive_La','relaxed_flip_La']}
ORDER6 = ('origin','adaptive_Ca','adaptive_La','amide_flip','relaxed_flip_Ca','relaxed_flip_La')
HA_TO_KCAL, EV_TO_KCAL = 627.509474, 23.06054783061903

def require(condition, message):
    if not condition: raise ValueError(message)

def pin(path):
    p = Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def check(record):
    require(pin(record['path']) == record, 'Changed pin: '+record['path'])
    return Path(record['path'])

def read(record): return json.loads((check(record) if isinstance(record,dict) else Path(record)).read_text())

def write(path, value):
    with Path(path).open('x') as f: json.dump(value,f,indent=2,allow_nan=False); f.write('\n')

def load(name, record):
    p = check(record); spec = importlib.util.spec_from_file_location(name,p)
    module = importlib.util.module_from_spec(spec); sys.modules[name] = module; spec.loader.exec_module(module)
    return module

def fixtures(m, base, old, corrected):
    inherited = base.parser_checks(m,old,corrected)
    archived = read(m['native_parent']); checks = {}
    keys = ('status','energy_hartree','printed_explicit_electrons','populations','group_indices',
            'group_charges','final_solver','final_trah_macro','scf_summary','checks')
    for record in m['mmol_collection_fixtures']:
        saved = read(record); tid, = saved['selected_task_ids']
        task, = [t for t in archived['tasks'] if t['task_id']==tid]
        pf = {'electronic_inventory':[base.electron_inventory(task,base.readxyz(check(task['xyz']),193))]}
        current = base.collect(archived,check(m['native_parent']),pf,old,corrected,{tid})
        before, = saved['cells']; after, = current['cells']
        check(before['output']); check(before['execution'])
        checks[tid] = all(before[k]==after[k] for k in keys) and all(after['checks'].values())
    require(len(checks)==6 and all(checks.values()), 'Six actual MMOL collection regressions failed')
    return {'passed':True,'inherited_ExaA_and_parser_checks':inherited,'six_actual_MMOL_checks':checks}

def sources(m, base):
    collection = read(m['source_collection']); preparation = read(m['source_preparation_collection'])
    row, = [r for r in collection['rows'] if r['source_case_id']==SOURCE]
    pr, = [r for r in preparation['rows'] if r['source_case_id']==SOURCE]
    source, = m['sources']; original = pr['original_row']
    require(collection['preparation']==m['source_preparation_collection'] and row['status']=='complete' and pr['status']=='prepared', 'Source not complete/prepared')
    require(source['source_case_id']==SOURCE and source['sample']==4 and (source['atom_count'],source['cap_count'])==(204,14), 'Wrong finite source')
    for k in ('case_id','case_key','record_id','selection_id','maps','state_signature','source_preparation','protonated_source'):
        require(source[k]==original[k], 'Original source lineage differs: '+k)
    for k, other in (('mapping','mapping'),('prepared_admission','admission'),('contact','contact'),('old_searches','old_searches')):
        require(source[k]==pr[other], 'Prepared lineage differs: '+k)
    require(source['new_proposals']==row['proposals'] and source['mapping']==row['prepared_mapping'] and source['prepared_admission']==row['admission'], 'Relaxed lineage differs')
    state=read(source['state_signature']); context=read(source['source_preparation']); admission=read(source['prepared_admission'])
    check(source['protonated_source']); check(source['mapping'])
    for key in ('origin','identities','source'): check(source['contact'][key])
    require(state['water_inventory']==[] and state['endpoint_charges']=={'Ca':-4,'La':-3} and state['multiplicities']=={'Ca':1,'La':1}, 'Physical chemical state differs')
    require(admission['pass'] and admission['status']=='prepared' and admission['caps_fixed'] and admission['all_physical_H_retained'] and admission['q0_paired_replay'], 'Rebased source admission failed')
    mapping=context['mapping']['source_to_qm']; caps=[x['qm_index'] for x in mapping if x['kind']=='sigma_link_H']
    require(sorted(x['qm_index'] for x in mapping)==list(range(1,204)) and len(caps)==len(state['caps'])==14, '204/14 physical mapping differs')
    require(read(source['maps']['Ca'])==read(source['maps']['La']), 'Paired original maps differ')
    coordinates={}
    for z in ('Ca','La'):
        for q in ORDER6:
            cell=row['matrix'][z][q]; require(cell['status']=='complete','Archived full pool unavailable')
            coordinates[z,q]=base.readxyz(check(cell['xyz']),204)
            require(coordinates[z,q][0][0]==z and [a[0] for a in coordinates[z,q]]==[a[0] for a in coordinates[z,'origin']], 'Candidate atom inventory/order changed')
            require(all(abs(float(coordinates[z,q][i][j])-float(coordinates[z,'origin'][i][j]))<=1e-12 for i in caps for j in (1,2,3)), 'Saved cap coordinates differ beyond inherited tolerance')
    for q in ORDER6:
        ca,la=coordinates['Ca',q],coordinates['La',q]
        require(ca[1:]==la[1:] and ca[0][1:]==la[0][1:], 'Archived paired candidate coordinates differ: '+q)
    composite=read(m['selection_collection']); sr,=[r for r in composite['rows'] if r['source_case_id']==SOURCE]
    for z in ('Ca','La'):
        for phase,pool in enumerate(('common3','common6')):
            expected=SELECTED[z][phase]
            for policy in ('operational','mathematical'):
                require(row['pools'][pool][policy]['selected_geometries'][z]==expected, 'Saved vacuum selection differs')
                picked=sr['pools'][pool][policy][z]
                require(picked['candidate']==expected and picked['components']['xyz']==row['matrix'][z][expected]['xyz'] and picked['components']['MACE']==row['matrix'][z][expected]['receipt'], 'Saved composite selection/XYZ differs')
    return row,pr,coordinates

def preflight(m,mp,base,old,corrected):
    import run_orca_task_manifest as runner
    import render_orca_runtime_input as runtime
    require(pin(runner.__file__)==m['implementation']['run_orca_task_manifest.py']['snapshot'], 'Wrong native executor import')
    require(pin(runtime.__file__)==m['implementation']['render_orca_runtime_input.py']['snapshot'], 'Wrong allocation renderer import')
    require(pin(sys.executable)==m['python'], 'Wrong Python executable')
    require(m['protocol_id']==PROTOCOL and m['expected_version']=='6.1.1' and m['La_ECP_core_electrons']==46 and m['selected_candidates']==SELECTED, 'Changed declared protocol')
    check(m['orca']); check(m['plan'])
    for record in m['runtime_files']: check(record)
    require(os.environ.get('METAL_ENV_WORKERS')=='1', 'One native worker required')
    policy=runtime.allocation_policy(workers=1,ranks=12)
    require(policy['workers']==1 and policy['ranks_per_worker']==12 and policy['allocated_cpus']==12 and policy['local_mpi_slots']==12, 'Actual single-node 12-MPI allocation required')
    job=subprocess.check_output(['scontrol','show','job',os.environ['SLURM_JOB_ID'],'-o'],text=True)
    fields=dict(v.split('=',1) for v in job.split() if '=' in v)
    require(fields.get('JobState')=='RUNNING' and fields.get('NumNodes')=='1' and fields.get('NumCPUs')=='12','Scheduler grant differs')
    runner.load_manifest_tasks(mp)
    expected={f'sample4_{z}_{q}_vacuum' for z in SELECTED for q in SELECTED[z]}
    require(len(m['tasks'])==4 and {t['task_id'] for t in m['tasks']}==expected,'Wrong four-cell denominator')
    row,prepared,coordinates=sources(m,base); inventory=[]; paths=[]; lineage=[]
    for t in m['tasks']:
        z,q=t['metal'],t['candidate']; cell=row['matrix'][z][q]; mr=read(t['saved_MACE']); request=read(t['saved_MACE_request'])
        require(t['source_xyz']==cell['xyz'] and t['xyz']['sha256']==cell['xyz']['sha256'] and t['saved_MACE']==cell['receipt'] and t['saved_forces']==cell['forces'],'Current saved cell binding differs')
        require(mr['request']==t['saved_MACE_request'] and request['xyz']['sha256']==cell['xyz']['sha256'],'Actual MACE request coordinates differ')
        check(request['xyz']); check(t['saved_forces']); check(request['manifest'])
        if 'manifest' in mr: require(mr['manifest']==request['manifest'],'MACE/request manifests differ')
        require(mr['status'] in ('computed','complete') and mr['energy_eV']==cell['energy_eV']==t['saved_MACE_eV'] and mr['forces']==cell['forces'],'Actual MACE energy/force receipt differs')
        require(request['charge']==t['charge'] and request['multiplicity']==t['multiplicity'],'Actual MACE request state differs')
        physical_state=mr.get('state',mr.get('input_state_check'))
        require(physical_state['atoms']==204 and physical_state['charge']==t['charge'] and physical_state['spin_multiplicity']==1,'Actual MACE state admission differs')
        prop=read(t['proposal']); point=prop['proposal']; new=t['phase']=='new'
        require(t['phase']==('new' if q.startswith('relaxed_flip_') else 'old'),'Selected phase label differs')
        require(request['manifest']==prop['manifest'] and request['task_id']==prop['task_id'],'Actual MACE request/proposal authorization differs')
        require(prop['status']=='proposal_available' and prop['metal']==z and prop['final_geometry']['physical_feasible'] and point['geometry_checks']['pass'],'Saved proposal admission failed')
        require(point['coordinate']['sha256']==cell['xyz']['sha256'] and point['MACE']==cell['receipt'] and point['forces']==cell['forces'],'Proposal point lineage differs')
        check(point['coordinate']); pm=read(prop['manifest']); task,=[v for v in pm['tasks'] if v['task_id']==prop['task_id']]
        source=m['sources'][0]
        require(task['charge']==t['charge'] and task['multiplicity']==1 and task['metal']==z and task['case_id']==prop['case_id'],'Proposal task state/source differs')
        if new:
            require(t['proposal']==row['proposals'][z] and task['mapping']==source['mapping'],'Rebased proposal mapping differs')
            contact=read(t['full_source_admission'])
            require(t['full_source_admission']==cell['full_source_admission'] and contact['pass'] and contact['task_id']==prop['task_id'] and contact['contact_inputs']==source['contact'] and contact['active_q_radian']==point['active_q_radian'],'Final full-source geometry admission differs')
        else:
            oldsearch=prepared['old_searches'][z]
            require(t['proposal']==oldsearch['result'] and prop['manifest']==oldsearch['manifest'] and task==oldsearch['task'] and task['mapping']==source['maps'][z],'Original proposal mapping differs')
            require(t['full_source_admission'] is None,'Invented old full-source admission')
        require(task['source_preparation']==source['source_preparation'],'Proposal source preparation differs')
        inp=check(t['input']); xyz=check(t['xyz']); count=base.readxyz(xyz,204)
        require(count==coordinates[z,q] and t['atom_count']==204 and t['sample']==4 and t['medium']=='vacuum' and t['charge']==(-4 if z=='Ca' else -3) and t['multiplicity']==1,'Native endpoint physical state differs')
        require(inp.read_text()==f"! r2SCAN-3c NoAutostart DefGrid3 TightSCF\n* xyzfile {t['charge']} 1 core.xyz\n",'Scientific input changed')
        require(xyz==inp.parent/'core.xyz' and Path(t['output_path'])==inp.parent/'endpoint.out','Task files are not isolated')
        paths.append(str(inp.parent)); inventory.append(base.electron_inventory(t,count)); lineage.append({'task_id':t['task_id'],'proposal':t['proposal'],'full_source_admission':t['full_source_admission'],'boundary_flag':prop['boundary_flag']})
    require(len(set(paths))==4,'Native directories overlap')
    mpi=shutil.which('mpirun'); ldd={name:subprocess.check_output(['ldd',str(Path(m['orca']['path']).parent/name)],text=True) for name in ('orca','orca_leanscf_mpi')}
    return {'status':'passed','manifest':pin(mp),'memory_policy':policy,'electronic_inventory':inventory,
        'source_checks':{'full_common6_paired_coordinates':True,'204_atoms_14_caps':True,'fixed_caps_tolerance_A':1e-12,'actual_old_and_new_selections':True,'lineage':lineage},
        'archived_parser_checks':fixtures(m,base,old,corrected),
        'runtime_observed':{'hostname':socket.gethostname(),'scontrol':job,'mpi':pin(mpi),'mpi_version':subprocess.check_output([mpi,'--version'],text=True),'affinity':sorted(os.sched_getaffinity(0)),
            'ldd':ldd,'libraries':[pin(p) for p in sorted({p for text in ldd.values() for p in re.findall(r'(/\S+)\s+\(0x',text)})],
            'environment':{k:v for k,v in os.environ.items() if k.startswith(('SLURM_','OMP','OMPI','MKL','OPENBLAS')) or k in ('PATH','LD_LIBRARY_PATH','LD_PRELOAD')}}}

def decision_gate(mp,path,sha):
    require(path is not None and sha,'Molecular execution requires the root post-transfer decision')
    record={'path':str(path.resolve()),'sha256':sha}; d=read(record)
    require(d.get('decision')=='execute_four_selected_Ca4_native_endpoints' and d.get('manifest')==pin(mp) and bool(d.get('interpretation')),'Missing positive exact-manifest decision')
    result=read(d['transfer_collection']); execution=read(d['transfer_execution'])
    require(result['protocol_id']=='uniform_terminal_amide_transfer_solvent_v1' and result['terminal'] and result['declared_contexts']==257 and result['declared_memberships']==260 and len(result['unique_context_rows'])==257 and len(result['rows'])==260,'Final uniform-transfer denominator unavailable')
    require(result['config']==execution['config'] and execution['status'] in ('complete','partial') and execution['action']=='execute','Transfer execution is not terminal')
    check(result['config'])
    return {'decision':record,'transfer_collection':d['transfer_collection'],'transfer_execution':d['transfer_execution'],'interpretation':d['interpretation']}

def scout_gate(m,mp,path,sha):
    require(path is not None and sha,'New-selected execution requires reviewed old-selected scouts')
    record={'path':str(path.resolve()),'sha256':sha}; d=read(record)
    require(d.get('decision')=='continue_two_new_selected' and d.get('manifest')==pin(mp) and bool(d.get('interpretation')),'Old-selected review differs')
    wanted={t['task_id'] for t in m['tasks'] if t['phase']=='old'}; observed=set()
    for record in d['scout_collections']:
        c=read(record); cell,=c['cells']
        require(c['manifest']==pin(mp) and c['complete_cells']==1 and cell['task_id'] in wanted and cell['status']=='complete_diagnostic' and all(cell['checks'].values()),'Old-selected scout not admitted')
        check(cell['output']); check(cell['execution']); observed.add(cell['task_id'])
    require(observed==wanted,'Both old-selected scouts required')
    return {'decision':pin(path),'scout_collections':d['scout_collections'],'interpretation':d['interpretation']}

def collect(m,mp,pf,base,old,corrected,selected):
    result=base.collect(m,mp,pf,old,corrected,selected)
    result.update(schema='A0AC_Ca4_selected_native_collection.v1',protocol_id=PROTOCOL,declared_cells=4,
        scientific_scope='selected old/new metal-specific paths; no common perturbation, pool minimum, solvent or classifier')
    cells={c['task_id']:c for c in result['cells']}; works={model:{} for model in ('DFT','MACE')}; discrepancies={}
    for z in SELECTED:
        tasks={t['phase']:t for t in m['tasks'] if t['metal']==z}
        works['MACE'][z]=(tasks['new']['saved_MACE_eV']-tasks['old']['saved_MACE_eV'])*EV_TO_KCAL
        oldcell=cells.get(tasks['old']['task_id']); newcell=cells.get(tasks['new']['task_id'])
        works['DFT'][z]=(newcell['energy_hartree']-oldcell['energy_hartree'])*HA_TO_KCAL if oldcell and newcell and all(c['status']=='complete_diagnostic' for c in (oldcell,newcell)) else None
        discrepancies[z]=None if works['DFT'][z] is None else works['DFT'][z]-works['MACE'][z]
    contrast={model:(v['Ca']-v['La'] if all(x is not None for x in v.values()) else None) for model,v in works.items()}
    result.update(selected_candidates=SELECTED,per_metal_old_to_new_work_kcal_mol=works,paired_selected_change_kcal_mol=contrast,
        DFT_minus_MACE_per_metal_kcal_mol=discrepancies,DFT_minus_MACE_paired_change_kcal_mol=None if contrast['DFT'] is None else contrast['DFT']-contrast['MACE'],
        interpretive_resolution_kcal_mol=1.0,resolution_is_acceptance_threshold=False,resolution_is_DFT_error_bar=False,
        units={'DFT':'hartree','MACE':'eV','work':'kcal/mol','HA_TO_KCAL':HA_TO_KCAL,'EV_TO_KCAL':EV_TO_KCAL})
    return result

def main():
    p=argparse.ArgumentParser(); p.add_argument('action',choices=('preflight','execute','collect')); p.add_argument('manifest',type=Path); p.add_argument('--manifest-sha',required=True)
    p.add_argument('--decision',type=Path); p.add_argument('--decision-sha'); p.add_argument('--scout-review',type=Path); p.add_argument('--scout-review-sha'); p.add_argument('tasks',nargs='*'); a=p.parse_intermixed_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'All scientific checks/execution/collection require Slurm')
    mp=check({'path':str(a.manifest.resolve()),'sha256':a.manifest_sha}); m=read(mp)
    for entry in m['implementation'].values(): check(entry['snapshot'])
    require(pin(__file__)==m['implementation']['REFERENCE.py']['snapshot'],'Wrong adapter snapshot')
    known={t['task_id'] for t in m['tasks']}; selected=set(a.tasks) if a.tasks else known
    require(len(a.tasks)==len(set(a.tasks)) and selected<=known and (a.action!='execute' or len(a.tasks)==1),'Execution requires one explicit task; no unknown/duplicate IDs')
    label=next(iter(selected)) if len(selected)==1 else 'all_four'
    attempt=mp.parent/'attempts'/f"{os.environ['SLURM_JOB_ID']}_{a.action}_{label}"; attempt.mkdir(parents=True,exist_ok=False)
    receipt={'status':'started','manifest':pin(mp),'action':a.action,'selected_task_ids':sorted(selected),'started_utc':dt.datetime.now(dt.timezone.utc).isoformat()}
    try:
        if a.action=='execute':
            write(attempt/'ROOT_DECISION.json',decision_gate(mp,a.decision,a.decision_sha))
            task,=[t for t in m['tasks'] if t['task_id'] in selected]
            if task['phase']=='new': write(attempt/'SCOUT_REVIEW.json',scout_gate(m,mp,a.scout_review,a.scout_review_sha))
        with ExitStack() as stack:
            for t in sorted((t for t in m['tasks'] if t['task_id'] in selected),key=lambda t:t['task_id']):
                lock=stack.enter_context((Path(t['output_path']).parent/'RUN.lock').open('a')); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            base=load('_A0AC_native_parent',m['implementation']['MMOL_REFERENCE.py']['snapshot']); old,corrected=base.modules(m)
            pf=preflight(m,mp,base,old,corrected); write(attempt/'PREFLIGHT.json',pf)
            if a.action=='execute':
                from run_orca_task_manifest import run_manifest
                try: run_manifest(mp,orca_path=Path(m['orca']['path']),workers=1,nprocs=12,selected_tasks=selected,expected_runner_sha256=m['execution_policy']['task_runner']['sha256'])
                finally:
                    result=collect(m,mp,pf,base,old,corrected,selected); write(attempt/'COLLECTION.json',result)
                require(result['complete_cells']==1,'Native endpoint failed diagnostic collection')
            elif a.action=='collect': write(attempt/'COLLECTION.json',collect(m,mp,pf,base,old,corrected,selected))
            receipt['status']='complete'
    except BaseException as exc:
        receipt.update(status='failed',error=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc()); raise
    finally:
        receipt['finished_utc']=dt.datetime.now(dt.timezone.utc).isoformat(); write(attempt/'EXECUTION.json',receipt)

if __name__=='__main__': main()
