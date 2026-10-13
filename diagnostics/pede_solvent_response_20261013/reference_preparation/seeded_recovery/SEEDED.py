#!/usr/bin/env python3
"""Slurm-only single-endpoint seed checks and unchanged native execution."""
import argparse
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import re
import traceback
from COLLECT_V2 import check,collect as base_collect,gradient_criterion,original,pin,write

INPUT='! r2SCAN-3c NoAutostart DefGrid3 TightSCF\n%scf\n Guess MORead\n MOInp "initial.gbw"\n MaxIter 125\nend\n* xyzfile -3 1 core.xyz\n'

def read_proof(text):
    actual='\n'.join(line for line in text.splitlines() if not re.match(r'\s*\|\s*\d+>',line))
    checks={'MOREAD_header':bool(re.search(r'^\s*INITIAL GUESS:\s*MOREAD\s*$',actual,re.M)),
            'actual_seed_filename':bool(re.search(r'^\s*Guess MOs are being read from file:\s*initial\.gbw\s*$',actual,re.M)),
            'same_basis_read':'Input basis set matches current basis set (good)' in actual,
            'no_atomic_fallback':'INITIAL GUESS: MODEL POTENTIAL' not in actual}
    return {'checks':checks,'passed':all(checks.values()),'runtime_lines':[s for s in actual.splitlines() if any(x in s for x in ('INITIAL GUESS:','Guess MOs are being read from file:','Input basis set','Input geometry','Occupation numbers'))]}

def validate(m,mp):
    for r in m['implementation'].values():check(r['snapshot'])
    for r in m['runtime_files']:check(r)
    check(m['root_plan']);pm=check(m['primary_manifest']);primary=json.loads(pm.read_text())
    source=next(t for t in primary['tasks'] if t['task_id']=='sample1_Ca_water')
    target=next(t for t in primary['tasks'] if t['task_id']=='sample1_Ca_vacuum')
    if source!=m['source_task'] or target!=m['primary_task']:raise ValueError('Original source/target changed')
    if len(m['tasks'])!=1:raise ValueError('Exactly one endpoint required')
    t=m['tasks'][0]
    if (t['sample'],t['metal'],t['medium'],t['charge'],t['multiplicity'])!=(1,'Ca','vacuum',-3,1):raise ValueError('Target state mismatch')
    if source['charge']!=-3 or source['multiplicity']!=1 or source['metal']!='Ca':raise ValueError('Source state mismatch')
    for rec in m['seed_source'].values():check(rec)
    c=json.loads(Path(m['seed_source']['collection']['path']).read_text())
    row=next(r for r in c['cells'] if r['task_id']=='sample1_Ca_water')
    receipt=json.loads(Path(m['seed_source']['execution']['path']).read_text())
    if c['manifest']!=m['primary_manifest'] or row['status']!='complete_diagnostic' or not all(row['checks'].values()):raise ValueError('Seed source unaccepted')
    if row['execution']!=m['seed_source']['execution'] or row['output']!=m['seed_source']['output'] or m['seed_source']['gbw'] not in row['retained_artifacts']:raise ValueError('Seed source provenance mismatch')
    if receipt['task_id']!=source['task_id'] or receipt['manifest']!=m['primary_manifest'] or receipt['orca_executable']!=m['orca']:raise ValueError('Source execution state/runtime mismatch')
    if not (receipt['returncode']==0 and receipt['normal_termination'] and receipt['scf_converged']):raise ValueError('Source not terminal success')
    for rec in receipt['artifacts'].values():check(rec)
    gbw=check(m['seed_source']['gbw']);seed=check(t['initial_gbw'])
    if gbw.stat().st_mtime>dt.datetime.fromisoformat(receipt['finished_at_utc']).timestamp():raise ValueError('Source GBW modified after completion')
    if pin(seed)['sha256']!=pin(gbw)['sha256'] or seed.name!='initial.gbw':raise ValueError('Wrong staged GBW')
    for q in (source,target,t):check(q['xyz']);check(q['input'])
    if len({q['xyz']['sha256'] for q in (source,target,t)})!=1:raise ValueError('Coordinates/order changed')
    if Path(target['input']['path']).read_text()!='! r2SCAN-3c NoAutostart DefGrid3 TightSCF\n* xyzfile -3 1 core.xyz\n':raise ValueError('Vacuum primary policy changed')
    if Path(source['input']['path']).read_text()!='! r2SCAN-3c NoAutostart DefGrid3 TightSCF CPCM(Water)\n* xyzfile -3 1 core.xyz\n':raise ValueError('Seed source method/basis changed')
    if Path(t['input']['path']).read_text()!=INPUT:raise ValueError('Recovery changed more than explicit initial guess/declared125bound')
    module=original(m);atoms=module.readxyz(t['xyz']['path']);zs={'H':1,'C':6,'N':7,'O':8,'S':16,'Ca':20}
    if atoms[0][0]!='Ca' or sum(zs[a[0]] for a in atoms)+3!=838:raise ValueError('192-atom singlet838-electron inventory mismatch')
    for s in m['sources']:
        for key in ('prepared_admission','source_preparation','state_signature'):check(s[key])
    qualified=json.loads(check(m['actual_parser_recovery']).read_text())
    if qualified['status']!='complete' or qualified['code']['COLLECT_V2.py']['sha256']!=pin(Path(__file__).with_name('COLLECT_V2.py'))['sha256']:raise ValueError('Corrected collector unqualified')
    from run_orca_task_manifest import load_manifest_tasks
    load_manifest_tasks(mp)
    return module,{'status':'passed','manifest':pin(mp),'initial_gbw':pin(seed),'source_gbw':pin(gbw),
                   'electronic_inventory':[{'task_id':t['task_id'],'explicit_electrons':838,'all_electrons':838}],
                   'atom_count':192,'charge':-3,'multiplicity':1,'scientific_change':'initial orbitals only; native MaxIter125 made explicit'}

def tests(m):
    raw=check(m['moread_fixture']).read_text();water=check(m['seed_source']['output']).read_text()
    proof=read_proof(raw)
    checks={'real_MOREAD_fixture_passes':proof['passed'],'real_atomic_guess_rejected':not read_proof(water)['passed'],
            'echo_only_rejected':not read_proof('\n'.join('| 1> '+s for s in proof['runtime_lines']))['passed'],
            'wrong_read_filename_rejected':not read_proof(raw.replace('file: initial.gbw','file: wrong.gbw'))['passed'],
            'missing_basis_match_rejected':not read_proof(raw.replace('Input basis set matches current basis set (good)',''))['passed'],
            'fallback_rejected':not read_proof(raw+'\nINITIAL GUESS: MODEL POTENTIAL\n')['passed']}
    return {'checks':checks,'all_passed':all(checks.values()),'fixtures':[m['moread_fixture'],m['seed_source']['output']],
            'new_molecular_calls':0,'negative_fixture_scope':'in-memory parser text only'}

def collect(module,m,mp,pf):
    result=base_collect(module,m,mp,pf);r=result['cells'][0]
    if 'checks' in r:
        actual=check(r['output']).read_text();r['orbital_read_proof']=read_proof(actual)
        r['checks']['explicit_seed_read']=r['orbital_read_proof']['passed']
        r['checks']['native_MaxIter125']=re.findall(r'Maximum # iterations\s+MaxIter\s+\.{3,}\s*(\d+)',actual)==['125']
        r['checks']['vacuum_target']='CPCM SOLVATION MODEL' not in actual and 'CPCM Dielectric' not in actual
        r['checks']['same_printed_state']=bool(re.search(r'Total Charge\s+Charge\s+\.{3,}\s+-3\b',actual) and re.search(r'Multiplicity\s+Mult\s+\.{3,}\s+1\b',actual))
        r['checks']['immutable_initial_gbw']=pin(m['tasks'][0]['initial_gbw']['path'])['sha256']==m['seed_source']['gbw']['sha256']
        r['status']='complete_diagnostic' if all(r['checks'].values()) else 'failed_collection'
    result.update(schema='pede_conditional_same_state_seeded_collection.v1',declared_cells=1,
                  complete_cells=int(r['status']=='complete_diagnostic'),primary_denominator_unchanged=True,
                  source_initialization=m['seed_source'],interpretation='same-state initialization recovery; no ground-state or stability qualification')
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=('dry-run','test','execute','collect'));p.add_argument('manifest',type=Path);p.add_argument('sha');p.add_argument('output',type=Path);p.add_argument('--test-receipt',type=Path);p.add_argument('--test-sha');a=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('All molecular parsing/checks require Slurm')
    mp=check({'path':str(a.manifest.resolve()),'sha256':a.sha});m=json.loads(mp.read_text())
    if m['implementation']['SEEDED.py']['snapshot']!=pin(__file__):raise ValueError('Wrong adapter')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    status={'status':'started','action':a.action,'manifest':pin(mp),'driver':pin(__file__),'job_id':os.environ['SLURM_JOB_ID'],'requested_molecular_endpoints':0}
    try:
        with (mp.parent/'RUN.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            module,pf=validate(m,mp);write(out/'PREFLIGHT.json',pf)
            if a.action=='test':
                result=tests(m);write(out/'TESTS.json',result)
                if not result['all_passed']:raise ValueError('Real-artifact seed read checks failed')
                status['tests']=pin(out/'TESTS.json')
            if a.action=='execute':
                if not a.test_receipt or not a.test_sha:raise ValueError('Actual passed Slurm test receipt required')
                tr=json.loads(check({'path':str(a.test_receipt.resolve()),'sha256':a.test_sha}).read_text())
                if tr['status']!='complete' or tr['action']!='test' or tr['manifest']!=pin(mp) or tr['driver']!=pin(__file__):raise ValueError('Wrong test receipt')
                if not json.loads(check(tr['tests']).read_text())['all_passed']:raise ValueError('Unpassed tests')
                pp=Path(m['primary_task']['output_path']).with_name('endpoint.out.execution.json')
                if not pp.exists():raise ValueError('Primary not terminal; conditional recovery remains prepared only')
                primary=json.loads(pp.read_text())
                if primary['returncode']==0 and primary['normal_termination'] and primary['scf_converged']:raise ValueError('Primary converged; no automatic alternate seed')
                if primary['task_id']!=m['primary_task_id'] or primary['manifest']!=m['primary_manifest']:raise ValueError('Wrong primary failure receipt')
                with (Path(m['primary_manifest']['path']).parent/'RUN.lock').open('a') as plock:
                    fcntl.flock(plock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                    from run_orca_task_manifest import run_manifest
                    from render_orca_runtime_input import allocation_policy
                    pf['memory_policy']=allocation_policy(workers=1,ranks=12);write(out/'EXECUTION_PREFLIGHT.json',pf)
                    status['primary_failure']=pin(pp);status['requested_molecular_endpoints']=1
                    try:run_manifest(mp,orca_path=Path(m['orca']['path']),workers=1,nprocs=12,selected_tasks={m['tasks'][0]['task_id']},expected_runner_sha256=m['execution_policy']['task_runner']['sha256'])
                    finally:write(out/'COLLECTION.json',collect(module,m,mp,pf))
                    if json.loads((out/'COLLECTION.json').read_text())['complete_cells']!=1:raise ValueError('Seeded endpoint not admitted')
            if a.action=='collect':write(out/'COLLECTION.json',collect(module,m,mp,pf))
            status['status']='complete'
    except BaseException as exc:
        status['status']='failed';status['error']=f'{type(exc).__name__}: {exc}';status['traceback']=traceback.format_exc();raise
    finally:write(out/'EXECUTION.json',status)

if __name__=='__main__':main()
