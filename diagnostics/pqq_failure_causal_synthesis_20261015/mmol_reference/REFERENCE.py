#!/usr/bin/env python3
"""Six saved MMOL193 vacuum cells; pinned native executor and TRAH collector."""
import argparse
from collections import Counter
from contextlib import ExitStack
import datetime as dt
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import sys
import traceback

PROTOCOL='MMOL1770_Ca0_saved_common3_native_r2SCAN3c_vacuum_v1'
ATOMIC_NUMBERS={'H':1,'C':6,'N':7,'O':8,'S':16,'Ca':20,'La':57}

def pin(path):
    path=Path(path).resolve();h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(path),'sha256':h.hexdigest()}

def check(r):
    if pin(r['path'])['sha256']!=r['sha256']:raise ValueError('SHA mismatch: '+r['path'])
    return Path(r['path'])

def read(r):return json.loads((check(r) if isinstance(r,dict) else Path(r)).read_text())
def write(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,r):
    p=check(r);s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m

def modules(m):
    for r in m['implementation'].values():check(r['snapshot'])
    old=load('_exaa_pinned_pede_reference',m['implementation']['OLD_REFERENCE.py']['snapshot'])
    corrected=load('_exaa_corrected_Trah',m['implementation']['COLLECT_V2.py']['snapshot'])
    return old,corrected

def electron_inventory(t,rows):
    full=sum(ATOMIC_NUMBERS[r[0]] for r in rows)-t['charge'];explicit=full-(46 if t['metal']=='La' else 0)
    if explicit%2:raise ValueError('Singlet electron parity differs')
    return {'task_id':t['task_id'],'composition':dict(Counter(r[0] for r in rows)),
            'all_electrons':full,'explicit_electrons':explicit}

def readxyz(path, count=193):
    lines=Path(path).read_text().splitlines();n=int(lines[0]);rows=[x.split() for x in lines[2:] if x.strip()]
    if n!=count or len(rows)!=count:raise ValueError('Wrong declared atom count')
    if any(len(r)!=4 or not all(math.isfinite(float(v)) for v in r[1:]) for r in rows):raise ValueError('Invalid XYZ')
    return rows

def parser_checks(m,old,corrected):
    exaa=load('_mmol_exaa_reference_fixture',m['implementation']['EXAA_REFERENCE.py']['snapshot'])
    inherited=exaa.parser_checks(m,old,corrected)
    counted=load('_mmol_counted_collector_fixture',m['implementation']['COUNTED_COLLECT.py']['snapshot'])
    archived=read(m['fixture_manifest']);checks={};details=[]
    for fixture in m['collection_fixtures']:
        prior=read(fixture);tid,=prior['selected_task_ids']
        task,=[t for t in archived['tasks'] if t['task_id']==tid]
        view={**archived,'tasks':[task]}
        pf={'electronic_inventory':[electron_inventory(task,readxyz(check(task['xyz']),192))]}
        current=counted.corrected_collect(old,corrected,view,check(m['fixture_manifest']),pf)['cells'][0]
        saved,=prior['cells'];check(saved['output']);check(saved['execution'])
        keys=('status','energy_hartree','printed_explicit_electrons','populations','group_indices','group_charges','final_solver','final_trah_macro','scf_summary')
        ok=all(current[k]==saved[k] for k in keys) and all(current['checks'].values()) and all(saved['checks'].values())
        checks[tid]=ok;details.append({'task_id':tid,'collection':fixture,'status':current['status'],'native_checks':current['checks']})
    if len(checks)!=8 or not all(checks.values()):raise ValueError('Eight real ExaA count-aware collection regressions failed')
    return {'passed':True,'inherited':inherited,'eight_real_ExaA_checks':checks,'details':details,'scope':'saved real artifacts only; no native calls'}

def preflight(m,mp,old,corrected):
    import render_orca_runtime_input as runtime
    import run_orca_task_manifest as runner
    check(m['orca']);inventory=read(m['source_collection']);check(m['plan'])
    for r in m['runtime_files']:check(r)
    for r in m['implementation'].values():check(r['snapshot'])
    if pin(runner.__file__)!=m['implementation']['run_orca_task_manifest.py']['snapshot']:raise ValueError('Wrong imported native runner')
    if pin(runtime.__file__)!=m['implementation']['render_orca_runtime_input.py']['snapshot']:raise ValueError('Wrong imported allocation renderer')
    if m['protocol_id']!=PROTOCOL or m['expected_version']!='6.1.1' or m['La_ECP_core_electrons']!=46:raise ValueError('Changed native protocol')
    if os.environ.get('METAL_ENV_WORKERS')!='1':raise ValueError('One allocation-aware native worker required')
    policy=runtime.allocation_policy(workers=1,ranks=12)
    if policy['workers']!=1 or policy['ranks_per_worker']!=12 or policy['allocated_cpus']!=12 or policy['local_mpi_slots']!=12:raise ValueError('Actual single-node12MPI allocation required')
    job=subprocess.check_output(['scontrol','show','job',os.environ['SLURM_JOB_ID'],'-o'],text=True)
    fields=dict(x.split('=',1) for x in job.split() if '=' in x)
    if fields.get('JobState')!='RUNNING' or fields.get('NumNodes')!='1' or fields.get('NumCPUs')!='12':raise ValueError('Actual scheduler allocation differs')
    runner.load_manifest_tasks(mp)
    expected={f'sample{s}_{z}_{c}_vacuum' for s in (0,) for z in ('Ca','La') for c in ('origin','adaptive_Ca','adaptive_La')}
    if len(m['tasks'])!=6 or {t['task_id'] for t in m['tasks']}!=expected:raise ValueError('Wrong finite6cell matrix')
    if len(m['sources'])!=1 or {s['sample'] for s in m['sources']}!={0}:raise ValueError('Wrong target-only source')
    source_rows={};states=[];origins={}
    for source in m['sources']:
        selected,=[r for r in inventory['rows'] if r['case_id']==source['case_id']]
        if selected['source_case_id']!='mmol_1770-pqq-la_model__conditioned_Ca__seed-1_sample-0' or source['source_case_id']!=selected['source_case_id'] or source['role']!='target':raise ValueError('Wrong source role/sample')
        for key in ('prepared_admission','source_preparation','state_signature','maps','native_proposals','proposal_manifest'):
            if source[key]!=selected[key]:raise ValueError('Source provenance differs: '+key)
        a=read(source['prepared_admission']);c=read(source['source_preparation']);state=read(source['state_signature'])
        if a['status']!='prepared' or not a['mapping_eligible'] or not a['paired_maps_equal'] or (a['atom_count'],a['cap_count'])!=(193,13):raise ValueError('Source admission failed')
        if state['water_inventory']!=[] or state['endpoint_charges']!={'Ca':-3,'La':-2} or state['multiplicities']!={'Ca':1,'La':1}:raise ValueError('Chemical state differs')
        if a['state_signature']!=source['state_signature'] or a['maps']!=source['maps'] or a['representations']['context']['preparation']!=source['source_preparation']:raise ValueError('Admission/source graph joins differ')
        mapping=c['mapping']['source_to_qm']
        if sorted(x['qm_index'] for x in mapping)!=list(range(1,193)) or sum(x['kind']=='sigma_link_H' for x in mapping)!=13:raise ValueError('Nonmetal/cap mapping differs')
        if read(source['maps']['Ca'])!=read(source['maps']['La']):raise ValueError('Paired physical maps differ')
        for metal in ('Ca','La'):
            proposal=read(source['native_proposals'][metal]);task=selected['proposal_metadata'][metal]['task']
            if proposal['metal']!=metal or proposal['case_id']!=source['case_id'] or not proposal['final_geometry']['physical_feasible'] or not proposal['proposal']['geometry_checks']['pass']:raise ValueError('Saved donor proposal is not admitted')
            if task['case_id']!=source['case_id'] or task['mapping']!=source['maps'][metal] or task['source_preparation']!=source['source_preparation']:raise ValueError('Proposal task/source map differs')
            if proposal['manifest']!=source['proposal_manifest']:raise ValueError('Proposal manifest differs')
            if proposal['proposal']['coordinate']!=selected['matrix'][metal]['adaptive_'+metal]['xyz']:raise ValueError('Own-metal saved proposal geometry differs')
        check(source['proposal_manifest'])
        for z in ('Ca','La'):
            check(source['native_proposals'][z])
            if source['endpoints'][z]!=selected['matrix'][z]['origin']['xyz']:raise ValueError('Origin source binding differs')
            origins[(source['sample'],z)]=readxyz(check(source['endpoints'][z]))
        states.append(state);source_rows[source['sample']]=selected
    if states[0]['core_microstates']!=source_rows[0]['state']['core_microstates']:raise ValueError('Core chemical states differ')
    records=[];coordinates={};directories=[]
    for t in m['tasks']:
        source=source_rows[t['sample']];cell=source['matrix'][t['metal']][t['candidate']]
        if t['source_xyz']!=cell['xyz'] or t['xyz']['sha256']!=t['source_xyz']['sha256']:raise ValueError('Exact saved candidate bytes differ')
        mr=read(t['saved_MACE']);request=read(t['saved_MACE_request'])
        if t['saved_MACE']!=cell['MACE'] or mr['request']!=t['saved_MACE_request'] or request['xyz']!=t['source_xyz']:raise ValueError('Saved MACE/source request differs')
        if mr['status']!='complete' or mr['energy_eV']!=cell['MACE_eV'] or t['saved_MACE_eV']!=cell['MACE_eV'] or request['charge']!=t['charge'] or request['multiplicity']!=t['multiplicity']:raise ValueError('Saved MACE state/energy differs')
        rows=readxyz(check(t['xyz']));check(t['source_xyz']);inp=check(t['input']);q=-3 if t['metal']=='Ca' else -2
        if t['medium']!='vacuum' or t['charge']!=q or t['multiplicity']!=1 or inp.read_text()!=f'! r2SCAN-3c NoAutostart DefGrid3 TightSCF\n* xyzfile {q} 1 core.xyz\n':raise ValueError('Scientific template/state changed')
        if rows[0][0]!=t['metal'] or [x[0] for x in rows]!=[x[0] for x in origins[(t['sample'],t['metal'])]]:raise ValueError('Candidate composition/order changed')
        if Path(t['output_path']).parent.resolve()!=inp.parent.resolve() or Path(t['xyz']['path']).resolve()!=inp.parent/'core.xyz':raise ValueError('Task files not isolated in one directory')
        directories.append(str(inp.parent.resolve()));coordinates[(t['sample'],t['candidate'],t['metal'])]=rows;records.append(electron_inventory(t,rows))
    if len(set(directories))!=6:raise ValueError('Task output directories overlap')
    for s in (0,):
        for c in ('origin','adaptive_Ca','adaptive_La'):
            ca=coordinates[s,c,'Ca'];la=coordinates[s,c,'La']
            if ca[1:]!=la[1:] or ca[0][1:]!=la[0][1:]:raise ValueError('Ca/La finite-motion coordinates differ')
    mpi=shutil.which('mpirun')
    observed={'hostname':socket.gethostname(),'mpi':pin(mpi),'mpi_version':subprocess.check_output([mpi,'--version'],text=True),
      'affinity':sorted(os.sched_getaffinity(0)),'scontrol':job,
      'environment':{k:v for k,v in os.environ.items() if k.startswith(('SLURM_','OMP','OMPI','MKL','OPENBLAS')) or k in ('PATH','LD_LIBRARY_PATH','LD_PRELOAD')}}
    observed['ldd']={name:subprocess.check_output(['ldd',str(Path(m['orca']['path']).parent/name)],text=True) for name in ('orca','orca_leanscf_mpi')}
    observed['libraries']=[pin(x) for x in sorted({p for text in observed['ldd'].values() for p in re.findall(r'(/\S+)\s+\(0x',text)})]
    return {'status':'passed','manifest':pin(mp),'memory_policy':policy,'electronic_inventory':records,
      'runtime_observed':observed,'archived_parser_checks':parser_checks(m,old,corrected),'scientific_scope':'fixed saved donor motions in vacuum; no new optimization'}

def collect(m,mp,pf,old,corrected,task_ids):
    selected=[t for t in m['tasks'] if t['task_id'] in task_ids]
    view={**m,'tasks':selected};counted=load('_mmol_counted_collector',m['implementation']['COUNTED_COLLECT.py']['snapshot']);result=counted.corrected_collect(old,corrected,view,mp,pf);tasks={t['task_id']:t for t in selected}
    for row in result['cells']:
        t=tasks[row['task_id']];row.update(candidate=t['candidate'],source_role='target')
        if 'output' in row:
            text=check(row['output']).read_text(errors='replace')
            row['electronic_state_warning_lines']=[x.strip() for x in text.splitlines() if re.search(r'warning|non.?Aufbau|HOMO.*LUMO|negative.*gap|small.*gap',x,re.I)]
        if 'checks' in row:
            r=read(row['execution']);extra={'native_executable_pin':r['orca_executable']['sha256']==m['orca']['sha256'],
              'native_runner_pin':r['task_runner']['sha256']==m['execution_policy']['task_runner']['sha256'],
              'runtime_renderer_pin':r['runtime_renderer']['sha256']==m['execution_policy']['runtime_renderer']['sha256'],
              'actual12MPI':r['parallelism']['nprocs']==12 and r['parallelism']['omp_threads_per_rank']==1,
              'native_auxiliary_basis':'def2-mTZVPP/J' in text,
              'native_singlet':bool(re.search(r'Multiplicity\s+Mult\s+\.+\s+1\b',text)),
              'population_partition':sorted(i for values in row['group_indices'].values() for i in values)==list(range(t['atom_count'])),
              'template_identity':r['artifacts']['template_input']['sha256']==pin(check(t['input']))['sha256'],
              'coordinate_identity':r['artifacts']['xyz']['sha256']==pin(check(t['xyz']))['sha256']}
            row['checks'].update(extra);row['status']='complete_diagnostic' if all(row['checks'].values()) else 'failed_collection'
        row['density_response_qualification']='unestablished; electronic warnings retained, no stability/repeat/gradient qualification'
    result.update(schema='mmol_saved_common3_dft_collection.v1',protocol_id=PROTOCOL,source_collection=m['source_collection'],
      declared_cells=6,selected_cells=len(selected),selected_task_ids=sorted(task_ids),complete_cells=sum(r['status']=='complete_diagnostic' for r in result['cells']),
      scientific_scope='same-coordinate vacuum finite-motion reference; no solvent study, classifier fit, force qualification or production promotion')
    return result

def review_gate(m,mp,path,sha):
    if path is None or not sha:raise ValueError('Proposal execution requires explicit reviewed origin scout receipt')
    decision=read({'path':str(path.resolve()),'sha256':sha})
    if decision.get('decision')!='continue_four_proposals' or decision.get('manifest')!=pin(mp) or not decision.get('interpretation'):raise ValueError('Missing positive origin review for this exact manifest')
    observed=set()
    for record in decision['scout_collections']:
        c=read(record)
        if c['manifest']!=pin(mp) or c['complete_cells']!=1:raise ValueError('Origin scout collection unavailable')
        cell,=c['cells'];tid=cell['task_id']
        if tid not in m['scout_tasks'] or cell['status']!='complete_diagnostic' or not all(cell['checks'].values()):raise ValueError('Invalid scout for continuation')
        check(cell['execution']);check(cell['output']);observed.add(tid)
    if observed!=set(m['scout_tasks']):raise ValueError('Both origin scouts required')
    return {'decision':pin(path),'actual_origin_task_ids':sorted(observed),'interpretation':decision['interpretation']}

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=('preflight','execute','collect'));p.add_argument('manifest',type=Path);p.add_argument('--manifest-sha',required=True);p.add_argument('--scout-review',type=Path);p.add_argument('--scout-review-sha');p.add_argument('tasks',nargs='*');a=p.parse_intermixed_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('All checks/execution/collection require actual Slurm')
    mp=check({'path':str(a.manifest.resolve()),'sha256':a.manifest_sha});m=read(mp)
    if pin(__file__)!=m['implementation']['REFERENCE.py']['snapshot']:raise ValueError('Wrong adapter snapshot')
    known={t['task_id'] for t in m['tasks']};chosen=set(a.tasks) if a.tasks else known
    if len(a.tasks)!=len(set(a.tasks)) or not chosen<=known or (a.action=='execute' and len(a.tasks)!=1):raise ValueError('Execute requires exactly one explicit frozen task; other actions accept unique selected tasks')
    name=next(iter(chosen)) if len(chosen)==1 else 'selected-'+hashlib.sha256(','.join(sorted(chosen)).encode()).hexdigest()[:12]
    attempt=mp.parent/'attempts'/f"{os.environ['SLURM_JOB_ID']}_{a.action}_{name}";attempt.mkdir(parents=True,exist_ok=False)
    status={'status':'started','manifest':pin(mp),'action':a.action,'tasks':sorted(chosen),'started_utc':dt.datetime.now(dt.timezone.utc).isoformat()}
    try:
        with ExitStack() as stack:
            for t in sorted((t for t in m['tasks'] if t['task_id'] in chosen),key=lambda t:t['task_id']):
                lock=stack.enter_context((Path(t['output_path']).parent/'RUN.lock').open('a'));fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            old,corrected=modules(m);pf=preflight(m,mp,old,corrected);write(attempt/'PREFLIGHT.json',pf)
            if a.action=='execute':
                from run_orca_task_manifest import run_manifest
                tid=next(iter(chosen))
                if tid not in m['scout_tasks']:write(attempt/'SCOUT_REVIEW.json',review_gate(m,mp,a.scout_review,a.scout_review_sha))
                try:run_manifest(mp,orca_path=Path(m['orca']['path']),workers=1,nprocs=12,selected_tasks={tid},expected_runner_sha256=m['execution_policy']['task_runner']['sha256'])
                finally:
                    result=collect(m,mp,pf,old,corrected,chosen);write(attempt/(tid+'_COLLECTION.json'),result);write(attempt/'COLLECTION.json',result)
                if result['complete_cells']!=1:raise RuntimeError('Native endpoint did not pass collection: '+tid)
            elif a.action=='collect':write(attempt/'COLLECTION.json',collect(m,mp,pf,old,corrected,chosen))
            status['status']='complete'
    except BaseException as e:
        status.update(status='failed',error=f'{type(e).__name__}: {e}',traceback=traceback.format_exc());raise
    finally:
        status['finished_utc']=dt.datetime.now(dt.timezone.utc).isoformat();write(attempt/'EXECUTION.json',status)

if __name__=='__main__':main()
