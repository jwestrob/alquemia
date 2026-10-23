#!/usr/bin/env python3
"""Separate full-panel FD checks from preserved archived-u/v regression failure."""
import argparse,importlib.util,json,math,sys,time
from pathlib import Path
from RUN import pin,read,verify,save,joins

def collect(recovery_pin,out):
 t=time.monotonic();cpu=time.process_time();out.mkdir(parents=True,exist_ok=False);r=read(recovery_pin);inputs=read(r['inputs']);panel,reused=joins(inputs)
 result={'status':'started','implementation':pin(__file__),'recovery':recovery_pin,'inputs':r['inputs'],'panel':r['panel'],'checks':{},'full_panel_admitted':False,'original_gradient_regression_passed':False,'same_archived_uv_diagnostic':None,'root_final_admission':None,'error':None,'new_native_calls':0}
 try:
  assert r['status']=='completed_33_previously_uncalled_force_states' and r['error'] is None and r['new_native_program_attempts']==r['new_native_program_invocations']==33 and set(r['rows'])==set(inputs['uncalled_ids'])
  assert r['reused_rows']==reused and set(r['rows']).isdisjoint(reused);maps={**reused,**r['rows']};assert set(maps)==set(panel['states']) and len(maps)==102
  oldcfg=read(inputs['parent_configuration']);rows={};counts={'direct_initializations':0,'mutual_operator_applications':0,'baseline_gradient_calls':0,'induced_solves':0}
  for id,p in maps.items():
   row=read(p);slot=panel['states'][id];assert row['state']==slot and row['returncode']==0 and row['native_program_attempts']==row['native_program_invocations']==1
   if id!=inputs['retained_failed_center']:assert row['status']=='native_force_state_complete' and row['error'] is None
   else:assert row['status']=='failed' and row['error']=="ValueError('available zero regression channel failed')"
   assert row['native_binary']==oldcfg['native_binary'] and row['kernel']==oldcfg['kernel'] and row['combination_source']==oldcfg['combiner']
   for key in ('baseline_state','supplied_uv','force_bridge_result'):verify(row[key])
   for inputpin in slot['inputs'].values():verify(inputpin)
   assert row['minimum_corrected_eigenvalue']>0 and math.isfinite(row['solve_residual']) and math.isfinite(row['energy_kcal_mol']) and all(math.isfinite(x) for atom in row['gradient_kcal_mol_A'] for x in atom)
   assert row['raw_native']['supplied_uv_unchanged'] and row['counts']=={'direct_initializations':1,'mutual_operator_applications':3*len(slot['geometry']['xyz_A']),'baseline_gradient_calls':1,'induced_solves':0}
   for key in counts:counts[key]+=row['counts'][key]
   if row['zero_regression'] and id!=inputs['retained_failed_center']:assert row['zero_regression']['passed']
   rows[id]=row
  source=verify(inputs['parent_collector']);sys.path.insert(0,str(source.parent));spec=importlib.util.spec_from_file_location('original_force_slope',source);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
  for id,row in rows.items():
   s=row['state']
   if s['sign']!=1:continue
   mid=id[:-4]+'minus';cid=s['case_id']+'__'+s['coefficient_mode']+'__center';check=old.slope(row['energy_kcal_mol'],rows[mid]['energy_kcal_mol'],s['step_A'],rows[cid]['gradient_kcal_mol_A'],panel['directions'][s['case_id']][s['direction']]);check.update(case=s['case_id'],profile=s['coefficient_mode'],direction=s['direction'],step_A=s['step_A'],plus=maps[id],minus=maps[mid],center=maps[cid]);result['checks'][id[:-6]]=check
  assert len(result['checks'])==48
  passed=all(c['passed'] for c in result['checks'].values());result.update(status='FD_panel_passed_original_regression_unresolved' if passed else 'FD_panel_failed_original_regression_unresolved',FD_checks_passed=passed,native_states=102,reused_scout_states=6,reused_prior_successes=62,reused_computed_failed_center=1,new_states=33,native_calls_by_kind=counts,original_gradient_regression={'result':maps[inputs['retained_failed_center']],'record':rows[inputs['retained_failed_center']]['zero_regression']},zero_regressions={id:x['zero_regression'] for id,x in rows.items() if x['zero_regression']},cost_sources={'scout':inputs['scout_result'],'partial_prior':inputs['partial_result'],'remaining33':recovery_pin},limits='All48FD gates remain5e-5+1e-6*abs(analytic). The original donor gradient regression stays failed. This collector cannot admit the full panel: exact archived-u/v diagnostic plus root review remain required.')
 except Exception as exc:result['status']='failed_saved_collection';result['error']=repr(exc)
 result['cost']={'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'threads':1};save(out/'RESULT.json',result);return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--recovery',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=collect(pin(a.recovery),a.out);print(json.dumps({'status':r['status'],'checks':len(r['checks']),'error':r['error']}));raise SystemExit(0 if r.get('FD_checks_passed') else 1)
