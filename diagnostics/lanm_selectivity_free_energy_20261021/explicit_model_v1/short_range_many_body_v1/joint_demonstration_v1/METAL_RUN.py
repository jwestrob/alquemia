#!/usr/bin/env python3
"""One root-bound Slurm campaign. No native/quantum executables or adaptive retries."""
import argparse,concurrent.futures,hashlib,importlib.util,json,math,multiprocessing,os,platform,resource,subprocess,sys,time,traceback
from pathlib import Path

def pin(path):
 p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def checked(p):
 if pin(p['path'])!=p:raise ValueError('changed input '+p['path'])
 return Path(p['path'])
def read(p):return json.loads(checked(p).read_text())
def load(name,p):
 path=checked(p);spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def clean(x):
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v) for v in x]
 if isinstance(x,float) and not math.isfinite(x):return 'nonfinite:'+repr(x)
 if hasattr(x,'tolist'):return clean(x.tolist())
 return x
def write(path,obj):
 path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(clean(obj),indent=2,allow_nan=False)+'\n');os.replace(tmp,path)
def threads():
 for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
  if os.environ.get(k)!='1':raise ValueError('required one-thread environment '+k)
 import torch
 torch.set_num_threads(1);torch.set_num_interop_threads(1)
def allocation():
 job=os.environ.get('SLURM_JOB_ID','')
 if not job.isdigit() or os.environ.get('SLURM_CPUS_PER_TASK')!='24':raise ValueError('exact 24-CPU Slurm allocation required')
 p=subprocess.run(['/usr/local/bin/scontrol','show','job','-o',job],capture_output=True,text=True,check=True)
 fields=dict(tok.split('=',1) for tok in p.stdout.split() if '=' in tok)
 for k,v in {'Partition':'standard','QOS':'standard_expand','Account':'banfield','NumCPUs':'24','NumNodes':'1','NumTasks':'1','Requeue':'0'}.items():
  if fields.get(k)!=v:raise ValueError('unexpected allocation '+k+': '+str(fields.get(k)))
 return {'job_id':job,'host':platform.node(),'pid':os.getpid(),'scontrol':p.stdout,'allocated_CPUs':24,'workers':16,'threads_per_worker':1,'scheduler_parent_cost':'pending terminal sacct; not inferred from child CPU sums'}

WORK=None
def initialize_worker(sources,states,bundle,settings,out):
 global WORK
 threads();K=load('joint_worker_kernel',sources['kernel']);C=load('joint_worker_core',sources['core'])
 WORK=(C,K,states,bundle,settings,out)
def worker(index,start):
 C,K,states,bundle,settings,out=WORK;directory=Path(out)/('start_%02d'%index);directory.mkdir(exist_ok=False)
 def checkpoint(r):write(directory/'RESULT.json',dict(r,worker={'pid':os.getpid(),'host':platform.node(),'threads':1}))
 result=C.solve_start(index,start,states,bundle,K,settings,checkpoint)
 checkpoint(result);return {'index':index,'result':pin(directory/'RESULT.json'),'status':result['status'],'score':result.get('assessment',{}).get('score')}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--config',required=True);parser.add_argument('--release',required=True);parser.add_argument('--out',required=True);a=parser.parse_args()
 out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
 result={'status':'started','error':None,'configuration':pin(a.config),'release':pin(a.release),'implementation':pin(__file__),
         'new_native_calls':0,'new_quantum_calls':0,'new_geometry_searches':0,'start_results':[],'selected_candidate':None}
 write(out/'RESULT.json',result);code=1
 try:
  cfg=json.loads(Path(a.config).read_text());release=json.loads(Path(a.release).read_text());src=cfg['sources']
  if release.get('status')!='released_joint_demonstration_v2' or release['configuration']!=pin(a.config) or release['implementation']!=pin(__file__) or release['output']!=str(out):raise ValueError('exact root release/configuration/output required')
  if cfg['output']!=str(out) or src['runner']!=pin(__file__):raise ValueError('configuration identity')
  for p in src.values():checked(p)
  if read(src['declaration'])['revision']['version']!=2:raise ValueError('prospective v2 declaration required')
  comparison=read(src['model_comparison'])
  if comparison['parent']!=src['declaration'] or cfg['role']!='metal15':raise ValueError('prospective metal15 comparison authority required')
  settings=read(src['settings']);result['allocation']=allocation();threads()
  import numpy as np,scipy,torch
  from scipy.stats import qmc
  result['versions']={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'torch':torch.__version__}
  K=load('joint_parent_kernel',src['kernel']);C=load('joint_parent_core',src['core']);D=load('joint_parent_data',src['data'])
  bundle=D.load(read(src['data_inputs']),read(src['endpoint_map']));states,symmetry=C.prepare(bundle,K)
  write(out/'SYMMETRY.json',symmetry);write(out/'DATA_PROVENANCE.json',{'fixture_pins':bundle['fixture_pins'],'provenance':bundle['provenance'],'units':bundle['units'],'method_channel_rules':bundle['method_channel_rules'],'operator_count':len(states),'all_data_consumed':True})
  result.update(operator_count=len(states),symmetry=pin(out/'SYMMETRY.json'),data_provenance=pin(out/'DATA_PROVENANCE.json'));write(out/'RESULT.json',result)
  zero=np.zeros(len(C.ORDER));pred,records=C.full_evaluate(states,bundle,K,zero);baseline={'coefficients':zero.tolist(),'predictions':records,'assessment':C.assess(pred,bundle)}
  if min(r['minimum_normalized_eigenvalue'] for r in records.values())<settings['normalized_SPD_floor']:raise ValueError('baseline fails declared SPD floor')
  write(out/'BASELINE.json',baseline);result['baseline']=pin(out/'BASELINE.json')
  check=C.gradient_check(states,bundle,settings);write(out/'GRADIENT_CHECK.json',check);result['gradient_check']=pin(out/'GRADIENT_CHECK.json');write(out/'RESULT.json',result)
  if check['status']!='passed':raise ValueError('declared analytic objective gradient check failed; no starts released')
  design=qmc.Sobol(d=4,scramble=True,seed=20261023).random_base2(4)[:15]*.2-.1
  starts=[np.zeros(4).tolist(),*design.tolist()];assert len(starts)==16
  write(out/'STARTS.json',{'generator':'scipy.stats.qmc.Sobol(d=4,scramble=True,seed=20261023).random_base2(4)[:15] *.2-.1, prepended zero','starts':starts,'settings':src['settings']});result['starts']=pin(out/'STARTS.json');write(out/'RESULT.json',result)
  with concurrent.futures.ProcessPoolExecutor(max_workers=16,mp_context=multiprocessing.get_context('spawn'),initializer=initialize_worker,initargs=(src,states,bundle,settings,str(out))) as pool:
   futures={pool.submit(worker,i,start):i for i,start in enumerate(starts)}
   for future in concurrent.futures.as_completed(futures):
    try:row=future.result()
    except Exception as exc:row={'index':futures[future],'status':'worker_exception','error':repr(exc),'traceback':traceback.format_exc()}
    result['start_results'].append(row);result['start_results'].sort(key=lambda r:r['index']);write(out/'RESULT.json',result)
  candidates=[];Hcount=hcount=Vcount=Jcount=0;workerCPU=0.
  for r in result['start_results']:
   if 'result' not in r:continue
   q=read(r['result']);workerCPU+=q.get('cost',{}).get('CPU_seconds',0.)
   # Intent checkpoints also survive abnormal process termination; these counts
   # describe returned solver records, not an assertion that failed starts did no work.
   Hcount+=int('H_solver' in q);hcount+=int('h_solver' in q);Vcount+=int('V_solver' in q);Jcount+=int('joint_solver' in q)
   if q['status']=='complete_stable_candidate':candidates.append((q['assessment']['score'],q['index'],q,r['result']))
  result['returned_solver_records']={'H_SLSQP':Hcount,'h_BVLS':hcount,'V_BVLS':Vcount,'joint_SLSQP':Jcount};result['completed_worker_CPU_seconds']=workerCPU
  if candidates:
   _,index,best,bestpin=min(candidates,key=lambda r:r[:2]);selected=best['candidates'][best['selected_stage']]
   result['selected_candidate']={'start_index':index,'stage':best['selected_stage'],'result':bestpin,'coefficients':selected['coefficients'],'score':selected['assessment']['score'],'criterion':'minimum declared five-family score among exact-box/all64-SPD feasible conditional or joint endpoints; exact ties earlier start then conditional','physical_accuracy_gate_passed':None}
   v=np.zeros(len(C.ORDER));v[C.VIDX]=np.array(selected['coefficients'])[C.VIDX];p,r=C.full_evaluate(states,bundle,K,v)
   write(out/'V_ONLY_ABLATION.json',{'coefficients':v.tolist(),'refitted':False,'predictions':r,'assessment':C.assess(p,bundle)})
   result['V_only_ablation']=pin(out/'V_ONLY_ABLATION.json');result['status']='completed';code=0
  else:result['status']='completed_no_stable_candidate'
  result['interpretation']=['Finite consumed-data nonconvex development campaign; no global capacity or physical accuracy proof.','All64 operators constrained, including report-only Dy and eight La endpoints.','All8 trimer cycles reported; only4 W9 cycles fitted.','Water total forces unavailable from frozen caches; no claim of force preservation.','No native, quantum, geometry search, new basis, adaptive starts, or automatic continuation.']
 except Exception as exc:result.update(status='failed',error=repr(exc),traceback=traceback.format_exc())
 finally:
  usage=resource.getrusage(resource.RUSAGE_SELF);children=resource.getrusage(resource.RUSAGE_CHILDREN)
  result['cost']={'controller_wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'controller_maxRSS_KiB':usage.ru_maxrss,'reaped_children_user_seconds':children.ru_utime,'reaped_children_system_seconds':children.ru_stime,'scheduler_parent_allocation_cost':'root adds terminal sacct; prior fixture production excluded'}
  write(out/'RESULT.json',result)
 return code
if __name__=='__main__':raise SystemExit(main())
