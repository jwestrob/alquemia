"""Native OMOL stage adapter using established worker and receipt acceptance."""
from __future__ import annotations
import argparse
import fcntl
import json
import os
from pathlib import Path
import subprocess
import time
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_hybrid import EV_TO_KCAL,accepted_attempt,check_atoms
import mace_omol as omol


def validate(manifest):
    m=read_json(manifest);verify(m['agreement']);verify(m['source_DFT_manifest'])
    if m['stage']!='second_shell_fixed' or len(m['tasks'])!=20 or m['model']!=omol.model(verify(m['software'])):
        raise InvalidArtifact('second-shell OMOL scope/model mismatch')
    for pin in m['implementation'].values():verify(pin)
    expected={(c,v,z) for c in ('1H4I','4MAE','1GLG','1F6S','6IP9') for v in ('core','expanded') for z in ('Ca','La')}
    if {(t['case_id'],t['variant'],t['metal']) for t in m['tasks']}!=expected:raise InvalidArtifact('task inventory changed')
    for t in m['tasks']:
        if check_atoms(xyz(verify(t['xyz'])),t['charge'])!=t['state']:raise InvalidArtifact('state changed')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('task cache identity changed')
    return {'status':'pass','tasks':20}


def execute(manifest):
    if not os.environ.get('SLURM_JOB_ID'):raise InvalidArtifact('allocation required')
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);sm=read_json(verify(m['software']))
    root=mp.parent/'mace_execution';root.mkdir(exist_ok=True)
    with (root/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for t in m['tasks']:
            td=root/t['task_id'];td.mkdir(exist_ok=True);previous=sorted(td.glob('attempt_*'))
            if any(accepted_attempt(p,t,mp) is not None for p in previous):continue
            if previous:raise InvalidArtifact('preserved partial attempt requires explicit unchanged retry')
            dest=td/'attempt_001';dest.mkdir()
            command=[str(verify(sm['python'])),str(verify(m['implementation']['second_shell_mace.py'])),
                     'worker','--manifest',str(mp),'--task-id',t['task_id'],'--output',str(dest)]
            write_new(dest/'started.json',{'manifest':record(mp),'task_id':t['task_id'],'command':command,
                'memory_mode':'native','slurm_job_id':os.environ['SLURM_JOB_ID'],'started_unix':time.time()})
            start=time.monotonic()
            with (dest/'worker.log').open('x') as log:
                child=subprocess.run(['/usr/bin/time','-v','-o',str(dest/'resources.txt'),*command],stdout=log,stderr=subprocess.STDOUT,check=False)
            receipt={'manifest':record(mp),'task_id':t['task_id'],'returncode':child.returncode,
                'wall_seconds':time.monotonic()-start,'memory_mode':'native','allocated_cpus':int(os.environ['SLURM_CPUS_PER_TASK']),
                'slurm_job_id':os.environ['SLURM_JOB_ID'],'log':record(dest/'worker.log'),'resource_usage':record(dest/'resources.txt'),
                'result':record(dest/'result.json') if (dest/'result.json').exists() else None}
            receipt['allocated_core_seconds']=receipt['wall_seconds']*receipt['allocated_cpus'];write_new(dest/'receipt.json',receipt)
            accepted=accepted_attempt(dest,t,mp)
            print(json.dumps({'task_id':t['task_id'],'status':'complete' if accepted else 'failed','wall_seconds':receipt['wall_seconds']}),flush=True)
            if accepted is None:return {'status':'partial_failure','task_id':t['task_id']}
    return {'status':'complete'}


def collect(manifest,output):
    m=read_json(manifest);rows=[]
    for t in m['tasks']:
        accepted=None
        for d in sorted((Path(manifest).parent/'mace_execution'/t['task_id']).glob('attempt_*')):
            accepted=accepted_attempt(d,t,manifest)
            if accepted:break
        rows.append({'task_id':t['task_id'],'case':t['case_id'],'variant':t['variant'],'metal':t['metal'],
                     'status':'complete' if accepted else 'unavailable','accepted':accepted})
    cases=[]
    for case in ('1H4I','4MAE','1GLG','1F6S','6IP9'):
        contrasts={}
        for variant in ('core','expanded'):
            pair={r['metal']:r for r in rows if r['case']==case and r['variant']==variant}
            if all(r['status']=='complete' for r in pair.values()):
                contrasts[variant]=(pair['Ca']['accepted']['energy_eV']-pair['La']['accepted']['energy_eV'])*EV_TO_KCAL
            else:contrasts[variant]=None
        cases.append({'case':case,'core_R_kcal_mol':contrasts['core'],'expanded_R_kcal_mol':contrasts['expanded'],
                      'delta_R_kcal_mol':contrasts['expanded']-contrasts['core'] if all(v is not None for v in contrasts.values()) else None})
    result={'manifest':record(manifest),'rows':rows,'cases':cases,'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete',
            'reference':None,'calibrated_decision':None,'solvent':'native OMOL vacuum; no CPCM equality claimed'}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    for name in ('validate','execute'):
        q=s.add_parser(name);q.add_argument('--manifest',required=True)
    q=s.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    q=s.add_parser('worker');q.add_argument('--manifest',required=True);q.add_argument('--task-id',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('op')
    if op=='worker':r=omol.worker(**a,memory_mode='native')
    else:r=globals()[op](**a)
    print(json.dumps({'status':r.get('status')}))


if __name__=='__main__':main()
