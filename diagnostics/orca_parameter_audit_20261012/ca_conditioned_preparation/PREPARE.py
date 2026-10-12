"""Finite source-graph preparation using unchanged existing envelope builders."""
from pathlib import Path
import os
THREAD_KEYS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
for k in THREAD_KEYS:os.environ[k]='1'
import argparse,copy,json,resource,shutil,sys,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from collections import Counter
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
import motion_envelope as envelope
from union_triple_transfer import mapping_equivalence
from affordable_common import read_json,record,verify,write_new,xyz
import numpy as np,scipy,gemmi
PROTOCOL='Nikasha_Ca_source_frozen_La_envelope4p3_preparation_v1'

def finalized(r,g,anchor):
 r.update(selection_id=g['selection_id'],roles=['Ca_conditioned_fixed_canonical_envelope_transfer'],triple_ids=[],state_anchor_source=g['members'][0])
 if r['status']=='prepared_awaiting_group_check':
  if anchor['status']!='prepared' or not anchor.get('mapping_eligible'):
   r.update(status='state_anchor_unavailable',reason='declared original anchor unavailable')
  elif anchor['state_key']!=r['state_key'] or anchor['protein_key']!=r['protein_key']:
   r.update(status='source_state_mismatch',reason='selected identities/H parents/caps/bonds/charges/protein differ from original group anchor')
  else:r['status']='prepared'
 return r

def task(row,g,config,out,anchor):
 start=time.monotonic();r=envelope.prepare_task(row,g,config,out)
 r=finalized(r,g,anchor);r['preparation_wall_seconds']=time.monotonic()-start
 write_new(Path(out)/'cases'/r['case_id']/'ADMISSION.json',r)
 return r

def compare_probe(new,old):
 assert new['status']==old['status']=='prepared'
 for k in ('state_key','protein_key','atom_count','cap_count','new_La_charge','paired_maps_equal'):
  assert new[k]==old[k],('replay metadata mismatch',new['case_id'],k)
 maximum=0.;exact=True
 for z in ('Ca','La'):
  a=new['representations']['context']['endpoints'][z];b=old['representations']['context']['endpoints'][z]
  assert (a['charge'],a['multiplicity'])==(b['charge'],b['multiplicity'])
  ar,br=xyz(verify(a['xyz'])),xyz(verify(b['xyz']));assert [r[0] for r in ar]==[r[0] for r in br]
  delta=float(np.max(np.abs(np.asarray([x[1:] for x in ar])-np.asarray([x[1:] for x in br]))))
  assert delta<=1e-12;maximum=max(maximum,delta);exact &= a['xyz']['sha256']==b['xyz']['sha256']
  a_map=read_json(verify(new['maps'][z]));b_map=read_json(verify(old['maps'][z]))
  # Copied original-cap vectors are coordinates, subject to the existing1e-12 A rule.
  for section in ('core','context'):
   for left,right in zip(a_map[section]['core_links'],b_map[section]['core_links']):
    if left[1]=='cap' and right[1]=='cap':
     d=float(np.max(np.abs(np.asarray(left[4:5])-np.asarray(right[4:5]))))
     assert d<=1e-12;maximum=max(maximum,d);left[4]=right[4]
     d=float(np.max(np.abs(np.asarray(left[5])-np.asarray(right[5]))))
     assert d<=1e-12;maximum=max(maximum,d);left[5]=right[5]
  assert a_map['checks']['pass'] and b_map['checks']['pass']
  for field in ('bond_max_error_A','jacobian_max_error_A_per_unit'):
   assert abs(a_map['checks'][field]-b_map['checks'][field])<=1e-12
   a_map['checks'][field]=b_map['checks'][field]
  maximum=max(maximum,mapping_equivalence(a_map,b_map))
 return {'case_id':new['case_id'],'status':'passed','maximum_coordinate_difference_A':maximum,'XYZ_hashes_identical':exact,'state_protein_and_discrete_mapping_identical':True}

def prepare(output,workers):
 assert workers in (1,2)
 out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);started=time.time();clock=time.monotonic()
 inv_path=ROOT/'diagnostics/orca_parameter_audit_20261012/ca_conditioned_inventory/INVENTORY.json';inv=read_json(inv_path)
 old_path=ROOT/'workspaces/motion_envelope_20260923/run_v2/PREPARATION.json';old=read_json(old_path)
 m=envelope.validate(verify(old['selection']));source_inputs=read_json(verify(m['source_inputs']))
 original=read_json(verify(source_inputs['source_preparation']));sources={r['case_id']:r for r in original['cases']}
 groups={r['selection_id']:r for r in m['groups']};oldrows={(r['case_id'],r['selection_id']):r for r in old['cases']}
 requests=sorted((r for r in inv['rows'] if r['required_work']=='expanded_preparation_then_proposals_and_components'),key=lambda r:r['case_id'])
 reused=sorted((r for r in inv['rows'] if r['required_work']=='reuse_saved_ownthree'),key=lambda r:r['case_id'])
 assert len(requests)==107 and len(reused)==3 and len(inv['rows'])==125
 for r in requests+reused:
  assert sources[r['case_id']]['status']=='prepared';g=groups[r['existing_target_selection_id']]
  assert g['union']==r['target_union'];anchor=oldrows[g['members'][0],g['selection_id']]
  assert anchor['status']=='prepared' and anchor['mapping_eligible']
 # Pin actual imported repository code, including unchanged legacy helpers.
 imported={name:record(mod.__file__) for name,mod in list(sys.modules.items()) if getattr(mod,'__file__',None) and str(ROOT/'scripts') in str(mod.__file__) and Path(mod.__file__).is_file()}
 runtime={'python':record(sys.executable),'version':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'gemmi':gemmi.__version__,
  'modules':{n:record(v.__file__) for n,v in [('numpy',np),('scipy',scipy),('gemmi',gemmi)]},'thread_environment':{k:os.environ[k] for k in THREAD_KEYS}}
 manifest={'protocol_id':PROTOCOL,'driver':record(__file__),'plan':record(Path(__file__).with_name('PLAN.md')),'inventory':record(inv_path),
  'old_preparation':record(old_path),'old_selection':old['selection'],'source_preparation':source_inputs['source_preparation'],'config':m['config'],
  'requests':requests,'reused':reused,'all125':inv['rows'],'groups':[groups[s] for s in sorted({r['existing_target_selection_id'] for r in requests+reused})],
  'implementation':m['implementation'],'imported_source_modules':imported,'runtime':runtime,'workers':workers,
  'scout_rule':'lexicographically_first_new_prepared_source','prospective_first_request':requests[0]['case_id'],
  'new_energies':0,'new_optimizations':0,'new_protonation':False,'new_waters':False,'production_changed':False}
 snapshot=out/'implementation';snapshot.mkdir()
 for name in ('PREPARE.py','PLAN.md'):shutil.copyfile(Path(__file__).with_name(name),snapshot/name)
 manifest['immutable_driver']=record(snapshot/'PREPARE.py');manifest['immutable_plan']=record(snapshot/'PLAN.md')
 write_new(out/'manifest.json',manifest)
 validation=[]
 for request in reused:
  cid,sid=request['case_id'],request['existing_target_selection_id'];g=groups[sid]
  r=task(sources[cid],g,m['config'],out/'validation'/sid,oldrows[g['members'][0],sid])
  validation.append(compare_probe(r,oldrows[cid,sid]))
 write_new(out/'REPLAY_CHECK.json',{'status':'passed','checks':validation,'molecular_calls':0})
 results={}
 with ProcessPoolExecutor(max_workers=workers) as ex:
  jobs={}
  for request in requests:
   cid,sid=request['case_id'],request['existing_target_selection_id'];g=groups[sid]
   jobs[ex.submit(task,sources[cid],g,m['config'],out/'prepared'/sid,oldrows[g['members'][0],sid])]=cid
  for future in as_completed(jobs):
   r=future.result();results[r['case_id']]=r
   print(json.dumps({'case_id':r['case_id'],'status':r['status'],'reason':r.get('reason')}),flush=True)
 rows=[]
 for i in inv['rows']:
  cid,sid=i['case_id'],i['existing_target_selection_id']
  if cid in results:r=results[cid]
  elif i['required_work']=='reuse_saved_ownthree':r=copy.deepcopy(oldrows[cid,sid]);r['prepared_artifact_reused']=True
  else:r={'case_id':cid,'selection_id':sid,'source':sources[cid]['source'],'status':'prior_source_preparation_excluded','reason':i['source_preparation_reason']}
  rows.append(r)
 ready=sorted([r for r in results.values() if r['status']=='prepared'],key=lambda r:r['case_id'])
 scout=ready[0] if ready else None
 usage=resource.getrusage(resource.RUSAGE_SELF);child=resource.getrusage(resource.RUSAGE_CHILDREN)
 receipt={'started_epoch':started,'finished_epoch':time.time(),'wall_seconds':time.monotonic()-clock,'workers':workers,'numerical_threads_per_worker':1,
  'CPU_seconds':usage.ru_utime+usage.ru_stime+child.ru_utime+child.ru_stime,'self_peak_RSS_KiB':usage.ru_maxrss,'child_peak_RSS_KiB':child.ru_maxrss,'scientific_energy_calls':0}
 summary={'protocol_id':PROTOCOL,'manifest':record(out/'manifest.json'),'replay_check':record(out/'REPLAY_CHECK.json'),'rows':rows,
  'counts':{'denominator':125,'new_preparation_attempts':107,'new_prepared':len(ready),'saved_preparations_reused':3,'prior_source_exclusions':15,
   'statuses':dict(Counter(r['status'] for r in rows))},'future_scout':{'case_id':scout['case_id'],'selection_id':scout['selection_id'],
     'endpoints':scout['representations']['context']['endpoints'],'maps':scout['maps'],'state_signature':scout['state_signature']} if scout else None,
  'receipt':receipt,'new_molecular_calls':0,'production_changed':False}
 write_new(out/'PREPARATION.json',summary);write_new(out/'EXECUTION.json',receipt)
 print(json.dumps({k:v for k,v in summary.items() if k!='rows'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--workers',type=int,required=True);a=p.parse_args();prepare(a.output,a.workers)
