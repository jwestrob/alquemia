#!/usr/bin/env python3
"""Run only 33 uncalled declared states using the immutable parent run_state."""
import argparse,concurrent.futures,hashlib,importlib.util,json,multiprocessing,os,platform,resource,socket,subprocess,sys,time
from pathlib import Path

def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
 assert pin(p['path'])==p,('changed pin',p);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def fp(cfg):return hashlib.sha256(json.dumps({k:v for k,v in cfg.items() if k!='root_release'},sort_keys=True,separators=(',',':')).encode()).hexdigest()
def allocation():
 jid=os.environ.get('SLURM_JOB_ID');assert jid and jid.isdigit(),'actual Slurm job required'
 raw=subprocess.check_output(['scontrol','show','job','-o',jid],text=True);fields=dict(x.split('=',1) for x in raw.split() if '=' in x)
 assert fields['JobId']==jid and fields['JobState']=='RUNNING' and fields['Account']=='banfield' and fields['QOS']=='standard_expand' and fields['Partition']=='standard'
 assert fields['UserId'].endswith('('+str(os.getuid())+')') and fields['NumNodes']=='1' and fields['NumCPUs']=='8' and fields['NumTasks']=='1' and fields['CPUs/Task']=='8'
 nodes=subprocess.check_output(['scontrol','show','hostnames',fields['NodeList']],text=True).split();assert len(nodes)==1 and nodes[0].split('.')[0]==socket.gethostname().split('.')[0]
 return {'job_id':jid,'host':socket.gethostname(),'fields':fields,'raw_scontrol':raw,'allocated_CPUs':8,'numerical_workers':4,'threads_per_worker':1,'controller_role':'metadata only'}
def load_parent(inputs):
 p=verify(inputs['parent_runner']);sys.path.insert(0,str(p.parent));spec=importlib.util.spec_from_file_location('native_force_parent_'+inputs['parent_runner']['sha256'][:12],p);m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m

def initialize(inputs):
 global PARENT,PCFG,PANEL,KERNEL,COMBINER
 for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
 PARENT=load_parent(inputs);PARENT.torch.set_num_threads(1);PARENT.torch.set_num_interop_threads(1);PCFG=read(inputs['parent_configuration']);PANEL=read(inputs['panel']);KERNEL=PARENT.module(PCFG['kernel']);COMBINER=PARENT.module(PCFG['combiner'])
def worker(id,folder):
 slot=PANEL['states'][id];r=PARENT.run_state(PCFG,PANEL,slot,KERNEL,COMBINER,Path(folder));return pin(Path(folder)/'RESULT.json')

def joins(inputs):
 parent_cfg=read(inputs['parent_configuration']);panel=read(inputs['panel']);scout=read(inputs['scout_result']);prior=read(inputs['partial_result'])
 assert parent_cfg['implementation']==inputs['parent_runner'] and parent_cfg['panel']==inputs['panel'] and parent_cfg['protocol']==inputs['parent_protocol'] and parent_cfg['metadata_helper']==inputs['parent_metadata_helper']
 for p in (inputs['parent_runner'],inputs['parent_protocol'],inputs['parent_metadata_helper']):verify(p)
 assert scout['status']=='native_force_scout_complete' and len(scout['rows'])==6 and prior['status']=='failed' and len(prior['rows'])==63 and prior['configuration']==inputs['parent_configuration']
 for key in ('panel','native_binary','kernel','combiner','build'):assert scout[key]==prior[key]==parent_cfg[key]
 assert set(scout['rows']).isdisjoint(prior['rows']);called={**scout['rows'],**prior['rows']};failed=[]
 for id,p in called.items():
  r=read(p);assert r['state']==panel['states'][id] and r['returncode']==0 and r['native_program_invocations']==r['native_program_attempts']==1
  if r['status']=='failed':
   assert id=='matched_dry_acetate_mono_central_La__zero__center' and r['error']=="ValueError('available zero regression channel failed')" and r['raw_native']['supplied_uv_unchanged'];assert not r['zero_regression']['gradient']['passed'] and all(x['passed'] for x in r['zero_regression']['available_channels'].values());failed.append(id)
  else:assert r['status']=='native_force_state_complete' and r['error'] is None
 assert failed==[inputs['retained_failed_center']]
 new=[id for id in panel['states'] if id not in called];assert new==inputs['uncalled_ids'] and len(new)==33 and all(panel['states'][id]['case_id']=='matched_dry_acetate_mono_central_La' for id in new)
 assert len(panel['states'])==102 and len(called)==69
 return panel,called

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();cfg=json.loads(a.config.read_text());inputs=read(cfg['inputs']);release=read(cfg['root_release'])
 assert cfg['implementation']==pin(__file__) and cfg['collector']==pin(Path(__file__).with_name('COLLECT.py')) and cfg['wrapper']==pin(Path(__file__).with_name('RUN.sbatch')) and cfg['workers']==4 and cfg['maximum_new_native_programs']==33
 assert release['status']=='approved_remaining33_native_force_states' and release['configuration_fingerprint']==fp(cfg) and release['output']==str(a.out.resolve()) and release['implementation']==cfg['implementation'] and release['inputs']==cfg['inputs']
 alloc=allocation();panel,reused=joins(inputs);a.out.mkdir(parents=True,exist_ok=False);t=time.monotonic();cpu=time.process_time()
 result={'status':'started','implementation':pin(__file__),'configuration':pin(a.config),'root_release':cfg['root_release'],'inputs':cfg['inputs'],'panel':inputs['panel'],'allocation':alloc,'reused_rows':reused,'retained_original_regression_failure':reused[inputs['retained_failed_center']],'rows':{},'submitted_ids':[],'pending_ids':[],'new_native_program_attempts':0,'new_native_program_invocations':0,'error':None,'full_panel_admitted':False};save(a.out/'RESULT.json',result)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn'),initializer=initialize,initargs=(inputs,)) as pool:
   ids=inputs['uncalled_ids']
   for start in range(0,33,4):
    wave=ids[start:start+4];assert not set(wave)&set(result['submitted_ids']);result['submitted_ids'].extend(wave);result['pending_ids']=wave.copy();save(a.out/'RESULT.json',result)
    futures={pool.submit(worker,id,str(a.out/'states'/id)):id for id in wave}
    for future in concurrent.futures.as_completed(futures):
     id=futures[future]
     try:
      p=future.result();r=read(p);result['rows'][id]=p;result['new_native_program_attempts']+=r['native_program_attempts'];result['new_native_program_invocations']+=r['native_program_invocations']
      if r['status']!='native_force_state_complete':result['error']={'state':id,'reason':r['error']}
     except Exception as exc:result['error']={'state':id,'worker_failure':repr(exc),'final_native_attempt_count_may_be_unavailable':True}
     result['pending_ids'].remove(id);save(a.out/'RESULT.json',result)
    if result['error'] is not None:break
 except Exception as exc:result['error']={'controller_failure':repr(exc),'final_native_attempt_count_may_be_unavailable':True}
 result['status']='completed_33_previously_uncalled_force_states' if result['error'] is None and len(result['rows'])==33 else 'failed_remaining33_native_force_states'
 result['cost']={'wall_seconds':time.monotonic()-t,'controller_CPU_seconds':time.process_time()-cpu,'worker_processes_total_CPU_seconds':resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime+resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime,'reported_state_native_CPU_seconds':sum(read(p)['cost']['native_child_CPU_seconds'] for p in result['rows'].values()),'host':platform.node(),'scheduler_accounting':None,'CPU_cost_warning':'Worker process totals may include descendant time; do not add reported native CPU to worker totals.'};save(a.out/'RESULT.json',result);print(json.dumps({'status':result['status'],'new_native_programs':result['new_native_program_invocations'],'error':result['error']}));return int(result['error'] is not None)
if __name__=='__main__':raise SystemExit(main())
