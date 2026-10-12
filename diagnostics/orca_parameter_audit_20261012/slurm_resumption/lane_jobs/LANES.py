#!/usr/bin/env python3
"""Two independent source lanes; unchanged MACE implementation and shared gate."""
import argparse,fcntl,hashlib,importlib.util,json,os,re,socket,subprocess,sys,time
from pathlib import Path

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')
TERMINAL={'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'}
def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):
 if isinstance(p,dict):
  assert pin(p['path'])==p;p=p['path']
 return json.loads(Path(p).read_text())
def write(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(plan):
 p=read(plan);assert pin(__file__)==p['wrapper']
 m=read(p['manifest']);assert m['implementation']==p['implementation']
 assert pin(m['implementation']['path'])==m['implementation']
 return p,m
def module(p):
 spec=importlib.util.spec_from_file_location('frozen_cpu_batch',p['implementation']['path'])
 b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b);return b
def states(jobs):
 found=subprocess.run(['sacct','-n','-P','-j',','.join(str(j) for j in jobs),'--format=JobIDRaw,State,ExitCode'],text=True,capture_output=True,check=True)
 rows={}
 for line in found.stdout.splitlines():
  parts=line.split('|')
  if parts[0] in {str(j) for j in jobs}:rows[parts[0]]={'state':parts[1].split()[0].rstrip('+'),'exit_code':parts[2]}
 return rows
def allocation(expected):
 job=os.environ.get('SLURM_JOB_ID');assert job and job.isdigit(),'actual Slurm allocation required'
 query=subprocess.run(['scontrol','show','job','-o',job],text=True,capture_output=True,check=True).stdout
 def field(name):
  value=re.search(r'(?:^|\s)'+name+r'=([^\s]+)',query);assert value,name;return value.group(1)
 assert field('JobState')=='RUNNING' and field('NumNodes')=='1'
 assert int(field('NumCPUs'))==expected and int(field('NumTasks'))==expected
 affinity=sorted(os.sched_getaffinity(0));assert len(affinity)>=expected
 assert int(os.environ['SLURM_CPUS_ON_NODE'])==expected
 return {'job_id':job,'host':socket.gethostname(),'scheduler_record':query,'affinity':affinity,'allocated_CPUs':expected,
         'environment':{k:v for k,v in os.environ.items() if k.startswith('SLURM_')},
         'memory':'actual scheduler default in record; no memory-sharing claim',
         'frozen_worker_placement_None_superseded_by_this_receipt':True}
def prepare(batch,drain,out):
 batch=batch.resolve();out=out.resolve();d=read(drain);manifest=pin(batch/'manifest.json');m=read(manifest)
 assert d['status']=='drained' and d['manifest']==manifest
 for k in ('live_owned_processes','stopped_owned_processes','incomplete_optimizer_task_ids'):assert d[k]==[]
 assert m['lanes']==2 and m['threads_per_model']==7 and m['maximum_active_numerical_threads']==14
 with (batch/'RUN.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  out.mkdir(parents=True,exist_ok=False)
  attempt=batch/'attempts'/('%03d'%(1+len(list((batch/'attempts').iterdir()))));attempt.mkdir()
  for lane in (0,1):(attempt/f'lane{lane}').mkdir()
  old=read(batch/'FINAL_COLLECTION.json')
  p={'protocol':'Slurm_independent_frozen_Ca_source_lanes_v1','wrapper':pin(__file__),'manifest':manifest,
     'implementation':m['implementation'],'drain':pin(drain),'attempt':str(attempt),'batch':str(batch),
     'submission_file':str(out/'SUBMISSION.json'),'already_terminal_sources':old['terminal_sources'],
     'lanes':{str(lane):[c['case_id'] for c in m['cases'] if c['lane']==lane] for lane in (0,1)},
     'resource_scope':'two natural disjoint source partitions,7CPU each,14 aggregate below60CPU/user; standard QOS/account unchanged',
     'scientific_method_changed':False,'qualification':'same four real controls and shared release; no new source before all pass'}
  write(out/'PLAN.json',p);print(json.dumps(p,indent=2))
def coordinate(plan):
 p,m=load(plan);out=Path(plan).parent;attempt=Path(p['attempt'])
 assert not os.environ.get('SLURM_JOB_ID'),'coordinator is metadata-only login client'
 with (Path(p['batch'])/'RUN.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  write(out/'COORDINATOR_READY.json',{'pid':os.getpid(),'host':socket.gethostname(),'plan':pin(plan),'manifest':p['manifest'],'scientific_calls':0})
  jobs=None;released=False;failure=None;terminal={}
  try:
   while not Path(p['submission_file']).exists():time.sleep(1)
   submission=read(p['submission_file']);assert submission['plan']==pin(plan)
   jobs=[str(submission['lane_job_ids'][str(lane)]) for lane in (0,1)];assert len(set(jobs))==2 and all(j.isdigit() for j in jobs)
   write(attempt/'RUN_STARTED.json',{'manifest':p['manifest'],'coordinator_pid':os.getpid(),'lane_job_ids':jobs,'maximum_active_numerical_threads':14})
   while True:
    terminal=states(jobs)
    checks=[attempt/f'lane{i}'/'QUALIFICATION.json' for i in (0,1)]
    if not released:
     if any(v['state'] in TERMINAL for v in terminal.values()):failure='lane_terminated_before_shared_qualification'
     elif all(q.exists() for q in checks):
      if all(read(q)['pass'] for q in checks):
       write(attempt/'QUALIFICATION_RELEASE.json',{'checks':[pin(q) for q in checks],'all_four_pass':True});released=True
      else:failure='qualification_failed'
     if failure and not (out/'QUALIFICATION_FAILURE.json').exists():
      write(out/'QUALIFICATION_FAILURE.json',{'reason':failure,'states':terminal,'barrier_released':False})
      # Only these two registered owned startup jobs; no source optimizer can have begun.
      for job in jobs:
       if terminal.get(job,{}).get('state') not in TERMINAL:subprocess.run(['scancel',job],check=True)
    if all(terminal.get(job,{}).get('state') in TERMINAL for job in jobs):break
    time.sleep(2)
  except BaseException as exc:
   failure=f'{type(exc).__name__}: {exc}'
   if jobs and not released:
    for job in jobs:
     if terminal.get(job,{}).get('state') not in TERMINAL:subprocess.run(['scancel',job],check=False)
   raise
  finally:
   write(out/'COORDINATOR_FINISH.json',{'plan':pin(plan),'released':released,'failure':failure,'states':terminal,'scientific_calls':0,'finished_epoch':time.time()})
def lane(plan,number):
 p,m=load(plan);out=Path(plan).parent;attempt=Path(p['attempt']);placement=allocation(7)
 assert number in (0,1) and all(os.environ.get(k)=='1' for k in THREADS)
 assert read(out/'COORDINATOR_READY.json')['plan']==pin(plan)
 write(attempt/f'lane{number}'/'SLURM_ALLOCATION.json',{'plan':pin(plan),'placement':placement})
 command=[m['optimizer_python'],'-u',p['implementation']['path'],'lane','--manifest',p['manifest']['path'],'--attempt',str(attempt),'--lane',str(number)]
 started=time.monotonic();child=subprocess.Popen(command,start_new_session=True)
 b=module(p);identity=b.identity(child.pid) or {'pid':child.pid,'start_ticks':None,'exited_before_identity':True}
 write(attempt/f'lane{number}'/'LANE_START.json',{**identity,'pgid':os.getpgid(child.pid),'command':command,'host':socket.gethostname(),'job_id':placement['job_id']})
 code=child.wait()
 write(attempt/f'lane{number}'/'SLURM_FINISH.json',{'returncode':code,'wall_seconds':time.monotonic()-started,'placement':placement})
 return code
def collect(plan):
 p,m=load(plan);out=Path(plan).parent;placement=allocation(1);submission=read(p['submission_file'])
 assert submission['plan']==pin(plan)
 jobs=[str(submission['lane_job_ids'][str(i)]) for i in (0,1)];terminal=states(jobs)
 assert all(terminal.get(job,{}).get('state') in TERMINAL for job in jobs),'afterany collector must follow both exact jobs'
 with (Path(p['batch'])/'RUN.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  assert (out/'COORDINATOR_FINISH.json').exists(),'metadata supervisor did not finish; inspect rather than invent completion'
  error=None if all(terminal[j]['state']=='COMPLETED' and terminal[j]['exit_code']=='0:0' for j in jobs) else 'one_or_more_lane_jobs_failed'
  result=module(p).collect(Path(p['manifest']['path']),Path(p['attempt']),error)
  write(out/'COLLECTOR_RESULT.json',{'plan':pin(plan),'placement':placement,'lane_states':terminal,
      'collection':pin(Path(p['batch'])/'FINAL_COLLECTION.json'),'terminal':result['terminal'],
      'terminal_sources':result['terminal_sources'],'native_matrix_complete':result['native_matrix_complete'],
      'native_matrix_unavailable':result['native_matrix_unavailable'],'failure':error})
 return 0 if result['terminal'] and error is None else 1
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);s=parser.add_subparsers(dest='op',required=True)
 q=s.add_parser('prepare')
 for k in ('batch','drain','output'):q.add_argument('--'+k,type=Path,required=True)
 for op in ('coordinate','lane','collect'):
  q=s.add_parser(op);q.add_argument('--plan',type=Path,required=True)
  if op=='lane':q.add_argument('--lane',type=int,choices=(0,1),required=True)
 a=parser.parse_args()
 if a.op=='prepare':prepare(a.batch,a.drain,a.output)
 elif a.op=='coordinate':coordinate(a.plan.resolve())
 elif a.op=='lane':sys.exit(lane(a.plan.resolve(),a.lane))
 else:sys.exit(collect(a.plan.resolve()))
