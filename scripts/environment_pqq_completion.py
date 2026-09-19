"""Finite local completion observer: no scientific execution or queue changes."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import resource
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--manifest',required=True)
    p.add_argument('--job-id',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--poll-seconds',type=int,default=60)
    a=p.parse_args()
    if not a.job_id.isdigit() or a.poll_seconds<30:raise ValueError('explicit Slurm job and nonbusy poll required')
    manifest=Path(a.manifest).resolve();out=Path(a.output).resolve()
    if out.exists():raise FileExistsError(out)
    sys.path.insert(0,str(manifest.parent/'implementation'))
    from affordable_common import InvalidArtifact,read_json,record,verify,write_new
    from hydration_square import endpoint
    from environment_pqq_dft import collect
    m=read_json(manifest);start=time.monotonic()
    for pin in m['implementation'].values():verify(pin)
    active={'PENDING','RUNNING','COMPLETING','CONFIGURING','SUSPENDED','RESIZING','REQUEUED','REQUEUE_FED','REQUEUE_HOLD'}
    while True:
        r=subprocess.run(['sacct','-j',a.job_id,'-n','-P','--format=JobIDRaw,State,ElapsedRaw,AllocCPUS,MaxRSS'],capture_output=True,text=True)
        lines=[line.split('|') for line in r.stdout.splitlines()]
        parent=next((line for line in lines if line[0]==a.job_id),None)
        if r.returncode==0 and parent and parent[1].split()[0].rstrip('+') not in active:
            break
        time.sleep(a.poll_seconds)
    statuses=[]
    for task in m['tasks']:
        try:
            e=endpoint(record(task['output_path']),record(task['output_path']+'.execution.json'),task['xyz'],task['input'])
            statuses.append({'task_id':task['task_id'],'status':'complete','energy_hartree':e['energy_hartree']})
        except (InvalidArtifact,OSError,ValueError,KeyError) as error:
            statuses.append({'task_id':task['task_id'],'status':'unavailable','reason':str(error),
                             'available_artifacts':[record(q) for q in (task['output_path'],task['output_path']+'.execution.json') if Path(q).exists()]})
    success=all(s['status']=='complete' for s in statuses)
    collection=None;collection_error=None
    if success:
        target=out.with_name('collected_'+a.job_id+'.json')
        try:
            collect(manifest,target);collection=record(target)
        except (InvalidArtifact,OSError,ValueError,KeyError) as error:
            collection_error=str(error)
    usage=resource.getrusage(resource.RUSAGE_SELF);children=resource.getrusage(resource.RUSAGE_CHILDREN)
    result={'job_id':a.job_id,'manifest':record(manifest),'observer_implementation':record(__file__),
            'scheduler_state':parent[1],'all_endpoints_complete':success,'endpoints':statuses,
            'collection':collection,'collection_error':collection_error,
            'allocated_wall_seconds':int(parent[2]),'allocated_cpus':int(parent[3]),
            'allocated_core_seconds':int(parent[2])*int(parent[3]),'allocated_GPU_seconds':0,
            'scheduler_accounting_rows':lines,'observer_wall_seconds':time.monotonic()-start,
            'observer_self_and_children_CPU_seconds':usage.ru_utime+usage.ru_stime+children.ru_utime+children.ru_stime,
            'new_scientific_calls':0,'changes_to_running_job':False}
    write_new(out,result);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
