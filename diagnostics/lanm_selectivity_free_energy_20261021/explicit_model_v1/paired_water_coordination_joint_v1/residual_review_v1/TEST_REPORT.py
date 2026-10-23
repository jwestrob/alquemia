#!/usr/bin/env python3
"""Small synthetic terminal records; no actual17676 outcomes or native calls."""
import copy,hashlib,importlib.util,json,math,os,sys,time
from pathlib import Path
import numpy as np
H=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('joint_residual_report_fixture',H/'REPORT.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def score(r):
 r=np.asarray(r);vals=r[:8]/math.sqrt(.25*.5/4)
 energy={'balanced_exchange':dict(zip(('water_minus','water_central','water_plus','aquo9_central'),vals[:4])),
  'individual_radial_work':dict(zip(('La_minus','La_plus','Dy_minus','Dy_plus'),vals[4:]))}
 energy['paired_radial_work']={side:energy['individual_radial_work']['Dy_'+side]-energy['individual_radial_work']['La_'+side] for side in ('minus','plus')}
 controls={g:{k:{'error_kcal_mol':float(v),'passed':bool(abs(v)<=1)} for k,v in rows.items()} for g,rows in energy.items()}
 per={kind:{q:{'rigid_RMS_kcal_mol_A':0.,'maximum_body_component_error_kcal_mol_A':0.,'passed':True} for q in m.GEOMETRIES} for kind in ('La','Dy','paired')}
 return {'objective':float(r@r),'weighted_residuals':r.tolist(),'energy_controls':controls,'per_fixture':per,'development_checks_passed':all(v['passed'] for rows in controls.values() for v in rows.values()),'blocks':{'energy':float(r@r)/.25,'La':0.,'Dy':0.,'paired':0.}}
def fixture(zero=False):
 base=np.zeros(413);base[0]=.3;center=np.full(6,.5);J=np.zeros((413,6))
 if not zero:J[:6]=np.eye(6)*8
 def row(v):return {'parameters':v.tolist(),'normalized_parameters':v.tolist(),'status':'complete','score':score(base+J@(v-center))}
 origin=row(center);records={'origin':origin};keys=[];widths=[]
 for j in range(6):
  for sign,name in ((-1,'minus'),(1,'plus')):
   v=center.copy();v[j]+=sign*.005;k=f'd{j}_{name}';records[k]=row(v);keys.append(k)
  a,b=(records[k] for k in keys[-2:]);widths.append(b['normalized_parameters'][j]-a['normalized_parameters'][j])
 matrix=np.column_stack([(np.array(records[keys[2*j+1]]['score']['weighted_residuals'])-np.array(records[keys[2*j]]['score']['weighted_residuals']))/widths[j] for j in range(6)])
 d,_,rank,sv=np.linalg.lstsq(matrix,-base,rcond=1e-6);direction=d.copy()
 if np.max(np.abs(direction))>.025:direction*=.025/np.max(np.abs(direction))
 trials=[];trialkeys=[];selected='origin'
 for name in ('scaled_GN','bounded_LS'):
  for frac,suffix in ((.5,'half'),(1.,'full')):
   v=center+frac*direction;k=name+'_'+suffix;records[k]=row(v);trialkeys.append(k);delta=v-center;pred=base+matrix@delta
   trials.append({'id':k,'parameters':v.tolist(),'normalized_delta':delta.tolist(),'predicted_objective':float(pred@pred)})
   if records[k]['score']['objective']<records[selected]['score']['objective']:selected=k
 entry={'status':'complete','outer_step':1,'origin_key':'origin','selected_key':selected,'derivative_keys':keys,'derivative_actual_widths':widths,
  'trial_keys':trialkeys,'trial_slots':trials,'Jacobian_normalized':matrix.tolist(),'linear_directions':{'rank':int(rank),'singular_values':sv.tolist(),'unrestricted_GN_not_evaluated':d.tolist()}}
 return {'status':'bounded_joint_coordination_water_complete','job_id':'17676','error':None,'candidate_frozen':True,'vectors':records,'initial':origin,'selected':records[selected],
  'family':{'parameters':[{'name':'p'+str(j),'lower':0.,'upper':1.} for j in range(6)]},'iterations':[entry],'completed_outer_steps':1,
  'stopping_reason':'no_actual_improvement_at_complete_trial_batch' if zero else 'all_joint_controls_passed','new_parameter_vectors':16,'new_native_program_invocations':0}
def main(out):
 t=time.monotonic();cpu=time.process_time();m.require(all(os.environ.get(k)=='1' for k in m.THREADS),'One thread required');out=Path(out);out.mkdir(parents=True,exist_ok=False);checks=[]
 def ok(v,label):
  if not v:raise AssertionError(label)
  checks.append(label)
 def rejects(r,label):
  try:m.summarize(r)
  except ValueError:checks.append(label)
  else:raise AssertionError(label)
 f=fixture();a=m.summarize(f);ok(a['selected']['energy_total']==10 and a['selected']['force_total']==15,'Exact control denominators');ok(a['selected']['energy_passed_count']==10 and a['selected']['force_passed_count']==15,'Synthetic actual pass retained')
 j=a['last_Jacobian'];ok(j['rank']==6 and not j['origin_equals_selected'] and not j['selected_point_projection_available'],'Different origin has no selected projection');ok(j['projection_applies_to']=='last_recorded_origin_only','Origin label explicit');ok(j['unrestricted_linear_projection']['weighted_loss']<1e-25,'Full-rank linear projection exact');ok(len(j['actual_trials'])==4 and a['steps'][0]['actual_improvement']>0,'Complete actual trial improvement retained')
 zero=m.summarize(fixture(True));j=zero['last_Jacobian'];ok(j['rank']==0 and j['condition'] is None and j['origin_equals_selected'],'Zero singular values handled with exact selected origin');ok(j['unrestricted_linear_projection']['weighted_loss']==j['origin_decomposition']['weighted_loss'],'Zero span does not invent projected improvement');ok(not j['global_family_impossibility_claimed'],'No global limitation claim')
 no=fixture(True);no['iterations']=[];no['completed_outer_steps']=0;no['selected']=no['initial'];a=m.summarize(no);ok(a['last_Jacobian']['status']=='unavailable','No-J controls summary retained')
 invalid=fixture();invalid['iterations'][0]['Jacobian_normalized'][0][0]+=1;rejects(invalid,'Reject mismatched stored Jacobian')
 invalid=fixture();invalid['selected']['score']['energy_controls']['balanced_exchange'].pop('water_minus');rejects(invalid,'Reject missing control denominator')
 invalid=fixture();invalid['status']='running';rejects(invalid,'Reject live/nonterminal result')
 invalid=fixture();invalid['iterations'][0]['status']='trials';rejects(invalid,'Reject incomplete terminal batch')
 result={'status':'passed','reporter':m.pin(H/'REPORT.py'),'test_source':m.pin(__file__),'specification':m.pin(H/'SPEC.md'),'checks':checks,'checks_passed':len(checks),'synthetic_only':True,'actual17676_outcomes_opened':False,'native_calls':0,'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'threads':1}
 (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','checks':len(checks),'result':m.pin(out/'RESULT.json')}))
if __name__=='__main__':main(sys.argv[1])
