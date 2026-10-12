#!/usr/bin/env python3
"""Slurm placement only; executes the unchanged, already-qualified CPU batch."""
from pathlib import Path
import argparse,datetime,hashlib,json,os,re,socket,subprocess,sys,time

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
THREAD_KEYS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')

def pin(path):
 p=Path(path).resolve()
 return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def read(path):return json.loads(Path(path).read_text())
def verify(p):
 assert pin(p['path'])==p,('changed artifact',p['path'])
 return Path(p['path'])
def write(path,data):
 with Path(path).open('x') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')

def prepared(manifest,drain):
 m=read(manifest)
 assert m['maximum_active_numerical_threads']==14 and m['threads_per_model']==7 and m['lanes']==2
 for k in ('implementation','scientific_scope','optimizer_executable'):verify(m[k])
 assert pin(m['optimizer_python'])==m['optimizer_executable']
 d=read(drain)
 assert d['status']=='drained' and d['manifest']==pin(manifest),'drain receipt must name this exact batch'
 for k in ('live_owned_processes','stopped_owned_processes','incomplete_optimizer_task_ids'):
  assert d[k]==[],('drain incomplete',k)
 return m,d

def placement():
 job=os.environ.get('SLURM_JOB_ID');assert job and job.isdigit(),'actual Slurm allocation required'
 command=['scontrol','show','job','-o',job]
 record=subprocess.run(command,text=True,capture_output=True,check=True)
 raw=record.stdout
 def field(name):
  match=re.search(r'(?:^|\s)'+re.escape(name)+r'=([^\s]+)',raw)
  assert match,('missing scheduler field',name)
  return match.group(1)
 assert field('Partition')=='standard-shared','this launcher is for the declared shared CPU partition'
 assert field('JobState')=='RUNNING' and field('NumNodes')=='1'
 assert field('NumCPUs')=='14','refusing a whole-node or insufficient CPU allocation'
 assert field('NumTasks')=='14','requested task slots changed'
 affinity=sorted(os.sched_getaffinity(0));assert len(affinity)>=14,'batch-shell affinity cannot supply the two seven-thread workers'
 assert int(os.environ['SLURM_CPUS_ON_NODE'])==14,'unexpected allocated CPU count'
 node=field('NodeList')
 nr=subprocess.run(['scontrol','show','node','-o',node],text=True,capture_output=True,check=True)
 return {'job_id':job,'hostname':socket.gethostname(),'job_query_command':command,'job_record':raw,
         'node_record':nr.stdout,'batch_shell_cpu_affinity':affinity,'allocated_CPUs':14,
         'SLURM_environment':{k:v for k,v in os.environ.items() if k.startswith('SLURM_')},
         'memory_policy':'no memory flags; actual scheduler memory retained in job_record; RAM sharing is not presumed',
         'placement_note':'This actual allocation receipt supersedes frozen worker scheduler_allocation=None placement metadata only.'}

def main(args):
 manifest=Path(args.manifest).resolve();drain=Path(args.drain).resolve()
 if args.dry_run:
  m=read(manifest);verify(m['implementation']);verify(m['optimizer_executable'])
  ready=False
  if drain.exists():prepared(manifest,drain);ready=True
  print(json.dumps({'status':'dry_run_no_scientific_calls','manifest':pin(manifest),'implementation':m['implementation'],
    'drain_present_and_accepted':ready,'execution_ready':False,'remaining_gate':'actual Slurm allocation/affinity checked inside job',
    'command':[m['optimizer_python'],'-u',m['implementation']['path'],'execute','--manifest',str(manifest)]},indent=2))
  return 0
 assert all(os.environ.get(k)=='1' for k in THREAD_KEYS),'set parent numerical thread environment before launch'
 out=Path(args.receipts).resolve();out.mkdir(parents=True,exist_ok=True)
 job=os.environ.get('SLURM_JOB_ID','no_allocation');destination=out/('job_'+job+'_START.json')
 started=time.monotonic();receipt={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
   'launcher':pin(__file__),'manifest':pin(manifest),'scientific_calls_by_launcher':0}
 try:
  m,d=prepared(manifest,drain);receipt['drain']=pin(drain);receipt['placement']=placement()
  receipt['original_implementation']=m['implementation']
  receipt['command']=[m['optimizer_python'],'-u',m['implementation']['path'],'execute','--manifest',str(manifest)]
  receipt['status']='preflight_pass';write(destination,receipt)
 except Exception as exc:
  receipt.update(status='refused_before_model_start',reason=str(exc));write(destination,receipt)
  raise
 child=subprocess.run(receipt['command'],cwd=ROOT)
 wall=time.monotonic()-started
 write(out/('job_'+job+'_FINISH.json'),{'start':pin(destination),'returncode':child.returncode,
   'wall_seconds_within_wrapper':wall,'allocated_CPU_slot_seconds_within_wrapper':14*wall,
   'utilization_measured':False,'final_collection':pin(manifest.parent/'FINAL_COLLECTION.json') if (manifest.parent/'FINAL_COLLECTION.json').exists() else None})
 return child.returncode

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('manifest','drain','receipts'):p.add_argument('--'+name,required=True)
 p.add_argument('--dry-run',action='store_true')
 sys.exit(main(p.parse_args()))
