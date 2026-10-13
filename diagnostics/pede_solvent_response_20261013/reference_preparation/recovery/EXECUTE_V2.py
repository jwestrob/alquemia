#!/usr/bin/env python3
"""External corrected collection around the unchanged frozen native runner."""
import argparse
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import traceback
import COLLECT_V2 as corrected

def main():
    p=argparse.ArgumentParser()
    p.add_argument('manifest',type=Path);p.add_argument('manifest_sha')
    p.add_argument('regression_receipt',type=Path);p.add_argument('regression_sha')
    p.add_argument('driver_sha');p.add_argument('output',type=Path);p.add_argument('tasks',nargs='+')
    a=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Actual Slurm required')
    corrected.check({'path':str(Path(__file__).resolve()),'sha256':a.driver_sha})
    mp=corrected.check({'path':str(a.manifest.resolve()),'sha256':a.manifest_sha})
    rp=corrected.check({'path':str(a.regression_receipt.resolve()),'sha256':a.regression_sha})
    rr=json.loads(rp.read_text());m=json.loads(mp.read_text())
    if rr['status']!='complete' or rr['new_molecular_calls']!=0:raise ValueError('Parser recovery is not successful collect-only evidence')
    regression=json.loads(corrected.check(rr['regression']).read_text())
    if not regression['all_passed'] or not all(regression['tests'].values()):raise ValueError('Regression checks did not pass')
    if corrected.pin(corrected.__file__)['sha256']!=rr['code']['COLLECT_V2.py']['sha256']:raise ValueError('Untested collector implementation')
    old_config=json.loads(corrected.check(rr['config']).read_text())
    original_manifest=json.loads(corrected.check(old_config['manifest']).read_text())
    if m['implementation']['REFERENCE.py']['snapshot']['sha256']!=original_manifest['implementation']['REFERENCE.py']['snapshot']['sha256']:
        raise ValueError('Frozen original collector/runner interface changed')
    if m['protocol_id']!=original_manifest['protocol_id']:raise ValueError('Scientific protocol changed')
    declared={t['task_id']:t for t in m['tasks']};frozen={t['task_id']:t for t in original_manifest['tasks']}
    if set(declared)!=set(frozen) or len(declared)!=8:raise ValueError('Wrong frozen matrix')
    for tid,t in declared.items():
        original_task=frozen[tid]
        for key in ('sample','metal','medium','charge','multiplicity'):
            if t[key]!=original_task[key]:raise ValueError('Task state changed: '+tid)
        for key in ('input','xyz','source_xyz'):
            if t[key]['sha256']!=original_task[key]['sha256']:raise ValueError('Scientific bytes changed: '+tid)
    if len(a.tasks)!=len(set(a.tasks)) or not set(a.tasks)<=set(declared):raise ValueError('Explicit unique frozen task IDs required')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    status={'status':'started','manifest':corrected.pin(mp),'parser_recovery':corrected.pin(rp),
            'driver':corrected.pin(__file__),'corrected_collector':corrected.pin(corrected.__file__),
            'job_id':os.environ['SLURM_JOB_ID'],'tasks':a.tasks,
            'started_utc':dt.datetime.now(dt.timezone.utc).isoformat()}
    module=None;pf=None
    try:
        with (mp.parent/'RUN.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            module=corrected.original(m);pf=module.preflight(m,mp)
            corrected.write(out/'PREFLIGHT.json',pf)
            from run_orca_task_manifest import run_manifest
            try:
                for tid in a.tasks:
                    run_manifest(mp,orca_path=Path(m['orca']['path']),workers=1,nprocs=12,selected_tasks={tid},
                                 expected_runner_sha256=m['execution_policy']['task_runner']['sha256'])
                    result=corrected.collect(module,m,mp,pf);corrected.write(out/(tid+'_COLLECTION.json'),result)
                    row=next(r for r in result['cells'] if r['task_id']==tid)
                    if row['status']!='complete_diagnostic':raise RuntimeError('Stop before next cell: '+tid+' '+row['status'])
            finally:
                corrected.write(out/'COLLECTION.json',corrected.collect(module,m,mp,pf))
            status['status']='complete'
    except BaseException as exc:
        status['status']='failed';status['error']=f'{type(exc).__name__}: {exc}';status['traceback']=traceback.format_exc();raise
    finally:
        status['finished_utc']=dt.datetime.now(dt.timezone.utc).isoformat();corrected.write(out/'EXECUTION.json',status)

if __name__=='__main__':main()
