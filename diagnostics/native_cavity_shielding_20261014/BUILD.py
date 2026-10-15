#!/usr/bin/env python3
"""Compile the thin diagnostic wrapper using the qualified native build helper."""
import hashlib,importlib.util,json,os,subprocess,sys,time
from pathlib import Path

def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
 if pin(p['path'])!=p:raise ValueError('Changed input: '+p['path'])
 return Path(p['path'])
mp=Path(sys.argv[1]).resolve()
if pin(mp)['sha256']!=sys.argv[2]:raise ValueError('Manifest changed')
m=json.loads(mp.read_text());verify(m['implementation']['BUILD.py'])
if pin(__file__)!=m['implementation']['BUILD.py']:raise ValueError('Wrong builder')
job=os.environ.get('SLURM_JOB_ID');assert job
allocation=subprocess.check_output(['scontrol','show','job','-o',job],text=True)
a=dict(v.split('=',1) for v in allocation.split() if '=' in v)
assert a['JobState']=='RUNNING' and a['NumCPUs']=='1' and a['NumNodes']=='1'
for p in m['build']['pins']+[m['implementation']['REPLAY.f90'],m['helper']]:verify(p)
s=importlib.util.spec_from_file_location('qualified_native_build',verify(m['helper']));n=importlib.util.module_from_spec(s);s.loader.exec_module(n)
out=mp.parent/'execution';out.mkdir(exist_ok=False)
env={k:v for k,v in os.environ.items() if k not in ('LD_PRELOAD','PYTHONPATH','CPATH','LIBRARY_PATH','CMAKE_PREFIX_PATH')}
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',BLIS_NUM_THREADS='1',LD_LIBRARY_PATH=m['build']['library_directory'])
t=time.monotonic();receipt={'status':'failed','manifest':pin(mp),'job_id':job,'allocation':allocation,'molecular_calls':0}
try:
 bp,b=n.compile_wrapper(m,mp,out,env);receipt.update(status='built',build=bp)
except BaseException as e:receipt['reason']=str(e);raise
finally:
 receipt['elapsed_seconds']=time.monotonic()-t
 (out/'EXECUTION.json').write_text(json.dumps(receipt,indent=2)+'\n')
