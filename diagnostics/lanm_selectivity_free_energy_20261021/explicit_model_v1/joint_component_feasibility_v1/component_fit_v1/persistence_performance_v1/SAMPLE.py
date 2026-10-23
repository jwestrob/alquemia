"""At most three nonblocking operational samples; no locals or fit-result reads."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

PID=3944417
BASE=Path(__file__).resolve().parent
OUT=BASE/'samples_v1'
PYSPY='/home/jwestrob/.local/bin/py-spy'
THREAD_KEYS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
for key in THREAD_KEYS:os.environ[key]='1'

def pin(path):
 path=Path(path).resolve();return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
def write(path,obj):
 with path.open('x') as f:json.dump(obj,f,indent=2,sort_keys=True);f.write('\n')
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def proc_stat(pid):
 p=Path('/proc')/str(pid)
 if p.stat().st_uid!=os.getuid():raise ValueError('Not own UID')
 raw=(p/'stat').read_text();end=raw.rfind(')');front=raw[:end+1];v=raw[end+2:].split()
 return {'pid':pid,'comm':front[front.index('(')+1:-1],'state':v[0],'ppid':int(v[1]),
         'user_ticks':int(v[11]),'system_ticks':int(v[12]),'threads':int(v[17]),'start_ticks':int(v[19]),
         'virtual_bytes':int(v[20]),'rss_pages':int(v[21])}
def snapshot():
 root=proc_stat(PID)
 io={}
 try:
  for line in (Path('/proc')/str(PID)/'io').read_text().splitlines():
   k,v=line.split(':',1);io[k]=int(v)
 except OSError as e:io={'error':repr(e)}
 processes={}
 for p in Path('/proc').iterdir():
  if p.name.isdigit():
   try:
    s=proc_stat(int(p.name));processes[s['pid']]=s
   except (OSError,ValueError,IndexError):pass
 descendants=[]
 for p,s in processes.items():
  if p==PID:continue
  parent=s['ppid'];seen=set()
  while parent and parent not in seen:
   if parent==PID:
    descendants.append({**s,'native_analyze_or_testgrad':s['comm'].startswith(('analyze','testgrad'))});break
   seen.add(parent)
   parent=processes.get(parent,{}).get('ppid',0)
 return {'utc':stamp(),'controller':root,'io':io,'descendants':descendants,'descendant_count':len(descendants),
         'native_descendant_count':sum(x['native_analyze_or_testgrad'] for x in descendants),
         'clock_ticks_per_second':os.sysconf('SC_CLK_TCK'),'page_size':os.sysconf('SC_PAGE_SIZE'),
         'snapshot_note':'Read-only non-atomic own-UID process snapshot; child counts are instantaneous.'}

def main():
 OUT.mkdir(exist_ok=False)
 declaration={'status':'declared_operational_sampling_only','created_utc':stamp(),'controller_pid':PID,
  'question':'Is whole-ledger JSON serialization a material current controller bottleneck?',
  'argv':[PYSPY,'dump','--pid',str(PID),'--nonblocking','--json'],'maximum_samples':3,'between_samples_seconds':2,
  'timeout_seconds_per_sample':8,'no_locals':True,'stop_on_first_pyspy_failure':True,'no_pause_or_runtime_mutation':True,
  'allowed_reads':['nonblocking Python frame locations','own-UID /proc stat and controller io','frozen source'],
  'forbidden_reads':['fit scientific values','partial result selection','Python locals'],'threads':1,
  'source':pin(Path(__file__)),'profiler':pin(PYSPY)}
 write(BASE/'DECLARATION.json',declaration)
 start,cpu=time.monotonic(),time.process_time();before_children=os.times();samples=[];identity=None
 for number in range(1,4):
  if number>1:time.sleep(2)
  sample={'number':number,'started_utc':stamp()}
  try:
   sample['before']=snapshot()
   if identity is None:identity=sample['before']['controller']['start_ticks']
   if sample['before']['controller']['start_ticks']!=identity:raise ValueError('Controller PID reused')
  except Exception as e:
   sample.update(status='proc_snapshot_failed',error=repr(e));samples.append(sample);write(OUT/f'sample_{number}.json',sample);break
  sample_start=time.monotonic()
  try:
   cp=subprocess.run(declaration['argv'],capture_output=True,timeout=8,check=False)
   (OUT/f'stack_{number}.json').write_bytes(cp.stdout)
   (OUT/f'stack_{number}.stderr').write_bytes(cp.stderr)
   sample.update(returncode=cp.returncode,stack=pin(OUT/f'stack_{number}.json'),stderr=pin(OUT/f'stack_{number}.stderr'))
   sample['status']='captured' if cp.returncode==0 else 'pyspy_failed_no_retry'
  except subprocess.TimeoutExpired as e:
   (OUT/f'stack_{number}.json').write_bytes(e.stdout or b'')
   (OUT/f'stack_{number}.stderr').write_bytes(e.stderr or b'')
   sample.update(status='pyspy_timeout_no_retry',error=repr(e),stack=pin(OUT/f'stack_{number}.json'),stderr=pin(OUT/f'stack_{number}.stderr'))
  sample['pyspy_wall_seconds']=time.monotonic()-sample_start
  try:sample['after']=snapshot()
  except Exception as e:sample['after']={'error':repr(e)}
  samples.append(sample);write(OUT/f'sample_{number}.json',sample)
  if sample['status']!='captured':break
 end_children=os.times()
 result={'status':'operational_samples_complete' if all(s['status']=='captured' for s in samples) else 'operational_sampling_stopped_on_failure',
  'declaration':pin(BASE/'DECLARATION.json'),'implementation':pin(Path(__file__)),'samples':samples,'sample_count':len(samples),
  'controller_start_ticks':identity,'wall_seconds':time.monotonic()-start,'sampler_cpu_seconds':time.process_time()-cpu,
  'child_cpu_seconds':end_children.children_user+end_children.children_system-before_children.children_user-before_children.children_system,
  'threads':1,'native_calls':0,'scientific_results_read':0,'fit_mutations':0}
 write(OUT/'RESULT.json',result)
 print(json.dumps({'status':result['status'],'result':pin(OUT/'RESULT.json'),'sample_count':len(samples),'wall_seconds':result['wall_seconds'],'sampler_cpu_seconds':result['sampler_cpu_seconds'],'child_cpu_seconds':result['child_cpu_seconds']}))
if __name__=='__main__':main()
