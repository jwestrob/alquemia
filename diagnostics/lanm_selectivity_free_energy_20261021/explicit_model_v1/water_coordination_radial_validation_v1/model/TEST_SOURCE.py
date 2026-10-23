#!/usr/bin/env python3
"""Source/geometry and rejected-candidate controls; never read actual17676 output."""
import copy,importlib.util,json,os,sys,time
from pathlib import Path
H=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('fresh_radial_model',H/'RUN.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def main(out):
 t=time.monotonic();cpu=time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False);checks=[]
 m.require(all(os.environ.get(k)=='1' for k in m.THREADS),'One thread required')
 def ok(v,label):
  if not v:raise AssertionError(label)
  checks.append(label)
 def rejects(fn,label):
  try:fn()
  except (ValueError,KeyError):checks.append(label)
  else:raise AssertionError(label)
 cfg=json.loads((H/'CONFIG_SOURCE.json').read_text());manifest=m.metadata(cfg)
 ok(cfg['candidate_result'] is None and cfg['candidate_review'] is None and cfg['root_release'] is None,'Candidate/review/release remain null')
 for gid,row in manifest['rows'].items():
  ok(row['geometry']['atoms']==28 and row['input']['atoms']==28,'Atom count '+gid)
  ok(len(row['input']['types'])==len(row['input']['bonds_zero_based'])==28,'Complete graph '+gid)
  ok(row['geometry']['radial_step_A']==(-.1 if '_minus_' in gid else .1),'Fixed radial step '+gid)
  for p in (row['geometry']['xyz'],row['input']['xyz']):m.verify(p);ok(True,'Geometry pin '+gid+' '+p['path'])
 for key,val in (('maximum_native_programs',9),('maximum_workers',5),('new_FD_energies',1)):
  c=copy.deepcopy(cfg);c[key]=val;rejects(lambda:m.metadata(c),'Reject changed limit '+key)
 c=copy.deepcopy(cfg);c['executables']['testgrad']['sha256']='0'*64;rejects(lambda:m.metadata(c),'Reject changed executable')
 rejects(lambda:m.candidate_gate(cfg),'Unbound candidate fails before result read')
 destination=out/'unbound_bind';destination.mkdir();rejects(lambda:m.bind(H/'CONFIG_SOURCE.json',destination),'Actual unbound bind fails closed')
 ok(not (destination/'BOUND.json').exists() and not (destination/'inputs').exists(),'Unbound bind created no chemistry inputs')
 # Mock only the result-reader interface. No running/terminal fit artifact is
 # opened; every case deliberately fails before data/model binding.
 rp={'path':str(m.WORK/'paired_water_coordination_joint_v1/run_v1/execution/RESULT.json'),'sha256':'synthetic_result'}
 vp={'path':str(out/'synthetic_review.json'),'sha256':'synthetic_review'}
 candidate={'status':'bounded_joint_coordination_water_complete','job_id':'17676','error':None,'candidate_frozen':True,'development_checks_passed':True,
  'config':cfg['fit_config'],'implementation':cfg['fit_implementation'],'selected':{'parameters':[0.]*6,'model':{g:{} for g in m.TRAINING_CELLS},'score':{'development_checks_passed':True}}}
 review={'status':'selected_ten_raw_review_passed','result':rp,'selected_parameters':[0.]*6,'selected_cell_ids':list(m.TRAINING_CELLS),
  'complete':10,'new_native_calls':0,'error':None,'all_selected_raw_replays_exact':True,'exact_selected_score_replayed':True,'all_development_controls_passed':True,'implementation':m.pin(Path(__file__))}
 c=copy.deepcopy(cfg);c.update(candidate_result=rp,candidate_review=vp);original=m.read;reads=[]
 for key,val in (('status','running'),('job_id','other'),('error','failed'),('candidate_frozen',False),('development_checks_passed',False)):
  altered=copy.deepcopy(candidate);altered[key]=val
  def fake(p):
   reads.append(p)
   if p==rp:return altered
   if p==vp:return review
   raise AssertionError('Unexpected actual result read')
  m.read=fake;rejects(lambda:m.candidate_gate(c),'Reject candidate '+key)
 for key,val in (('status','partial'),('complete',9),('new_native_calls',1),('error','failure'),('selected_cell_ids',list(m.TRAINING_CELLS)[:-1]),('all_selected_raw_replays_exact',False),('exact_selected_score_replayed',False),('all_development_controls_passed',False),('selected_parameters',[1.]*6)):
  altered=copy.deepcopy(review);altered[key]=val
  def fake(p):
   reads.append(p)
   if p==rp:return candidate
   if p==vp:return altered
   raise AssertionError('Unexpected actual result read')
  m.read=fake;rejects(lambda:m.candidate_gate(c),'Reject review '+key)
 m.read=original;ok(all(p in (rp,vp) for p in reads),'Only synthetic result/review objects accessed')
 ok(not any(out.rglob('analyze.out')) and not any(out.rglob('testgrad.out')),'No native output or molecular calls')
 result={'status':'passed','implementation':m.pin(__file__),'model_source':m.pin(H/'RUN.py'),'config':m.pin(H/'CONFIG_SOURCE.json'),
  'geometry_manifest':cfg['geometry_manifest'],'checks':checks,'checks_passed':len(checks),'new_native_calls':0,'candidate_bound':False,
  'fit_outcomes_opened':False,'validation_outcomes_opened':False,'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'threads':1}
 m.save(out/'RESULT.json',result);print(json.dumps({'status':'passed','checks':len(checks),'result':m.pin(out/'RESULT.json')}))
if __name__=='__main__':main(sys.argv[1])
