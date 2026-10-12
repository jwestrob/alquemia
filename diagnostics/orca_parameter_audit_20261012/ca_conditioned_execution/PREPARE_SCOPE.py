"""Pin finite Ca transfer scope/bands; no numerical-library or model calls."""
from pathlib import Path
import argparse,collections,hashlib,json
ROOT=Path(__file__).resolve().parents[3]
def pin(path):
 p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(path):return json.loads(Path(path).read_text())
def checked(p):assert pin(p['path'])==p;return read(p['path'])
def write(p,x):
 with Path(p).open('x') as f:json.dump(x,f,indent=2);f.write('\n')
def prepare(config,output):
 cfg=read(config);paths={k:ROOT/v for k,v in cfg.items()};d={k:read(v) for k,v in paths.items()}
 prep=d['preparation'];orig=d['origin_inputs'];ref=d['reference'];prior=d['prior_collection'];contract=d['proposal_contract']
 assert orig['preparation']==pin(paths['preparation']) and len(orig['tasks'])==212
 assert ref['manifest']==prior['manifest'] and ref['calibration']==prior['calibration']
 assert all(x['status']=='available' and x['available']==25 for x in ref['calibration'].values())
 assert contract['model']==orig['model'] and d['CPU_qualification']['status']=='passed_scoped_CPU_force_and_two_proposal_qualification'
 rows=[];new=[];reused=[];excluded=[];old={r['source_case_id']:r for r in prior['rows'] if r['source_conditioning_metal']=='Ca'}
 tasks_by_case=collections.defaultdict(dict)
 for t in orig['tasks']:
  assert t['metal'] not in tasks_by_case[t['source_case_id']]
  tasks_by_case[t['source_case_id']][t['metal']]=t
  for k in ('xyz','mapping','state_signature','context_preparation','prepared_admission'):assert pin(t[k]['path'])==t[k]
 for p in prep['rows']:
  cid=p['case_id'];source=p['source'];row={'source_case_id':cid,'selection_id':p['selection_id'],
   'protein':source['root_case_id'],'biological_group':source['biological_group'],'expected_class':source['expected_class'],
   'source_conditioning_metal':'Ca','preparation_status':p['status'],'preparation_reason':p.get('reason')}
  assert source['source_conditioning_metal']=='Ca' and row['expected_class'] in ('Ca','La')
  if p['status']!='prepared':row['execution_status']='excluded';excluded.append(cid)
  elif p.get('prepared_artifact_reused'):
   r=old[cid];assert r['selection_id']==p['selection_id'] and r['pool']['status']=='available'
   for cell in r['components'].values():
    assert cell['accepted']
    for point in cell['inputs'].values():checked(point)
   row.update(execution_status='reuse_actual_corrected_pool',case_id=r['case_id'],prior_collection=pin(paths['prior_collection']),prior_row=r,
      state_signature=p['state_signature'],maps=p['maps'],source_preparation=p['representations']['context']['preparation'])
   reused.append(cid)
  else:
   ts=tasks_by_case[cid];assert set(ts)=={'Ca','La'} and ts['Ca']['case_id']==ts['La']['case_id']
   assert ts['La']['charge']-ts['Ca']['charge']==1 and ts['Ca']['multiplicity']==ts['La']['multiplicity']==1
   row.update(execution_status='new_native_pool_required',case_id=ts['Ca']['case_id'],tasks=ts,state_signature=p['state_signature'],
      maps=p['maps'],source_preparation=p['representations']['context']['preparation'])
   new.append(cid)
  rows.append(row)
 assert len(rows)==125 and len(new)==106 and len(reused)==3 and len(excluded)==16
 r=checked(ref['manifest']);sc=d['CPU_qualification']['cost'];proposals=d['CPU_qualification']['proposals']
 calls_per_source=sum(p['fresh_force_calls'] for p in proposals)+4
 # This is a one-202-atom-scout extrapolation at8threads; planned7threads is unmeasured.
 forwards=106*calls_per_source+4;seconds_per_call=sc['per_call_seconds']['sum']/sc['actual_new_forward_calls']
 costs={'status':'planning_extrapolation_not_benchmark_or_budget','basis':pin(paths['CPU_qualification']),
  'measured_scout_atoms':202,'measured_model_threads':8,'planned_model_threads':7,'planned_concurrent_model_workers':2,
  'measured_scout_calls':sc['actual_new_forward_calls'],'measured_scout_total_CPU_seconds':sc['total_CPU_seconds'],
  'assumed_new_source_calls_for_estimate':calls_per_source,'projected_calls_including_four_initial_qualification':forwards,
  'projected_sum_model_wall_seconds':forwards*seconds_per_call,'ideal_two_lane_model_hours_at_observed8thread_rate':forwards*seconds_per_call/2/3600,
  'limits':['only one202atom source qualified','new sources169–215atoms','7thread and simultaneousworker throughput not yet measured','optimizer trialcounts variable; maxiter200 is unchangedmethod not newbudget','all actual starts/calls/failures/timing must be collected'],
  'GFN_stage_inputs_per_backend_maximum':1272,'GFN_new_molecules_not_launched':True}
 scope={'protocol_id':'Nikasha_corrected_expanded_Ca_conditioned_transfer_scope_v1','driver':pin(__file__),'config':pin(config),
  'plan':pin(Path(__file__).with_name('PLAN.md')),'inputs':{k:pin(v) for k,v in paths.items()},'model':contract['model'],
  'settings':{'optimizer':contract['optimizer_settings'],'pool':contract['pool_settings'],'scoring':r['settings']},
  'reference':{'source':pin(paths['reference']),'frozen_calibration':ref['calibration'],'new_fit':False},
  'rows':rows,'counts':{'declared_Ca_sources':125,'biological_groups':25,'new_sources':106,'reused_sources':3,'prior_PQQ_exclusions':15,'new_preparation_exclusions':1,
   'new_origin_calls':212,'optimizer_starts':212,'opposite_metal_MACE_cells_maximum':212,'GFN_inputs_per_backend_maximum':1272},
  'execution_resources':{'model_workers':2,'threads_per_model_worker':7,'optimizer_parent_threads':1,'parent_numerics_only_while_own_child_idle':True,'aggregate_active_numerical_threads':14,'root_headroom_threads':2,'user_total_threads':16},
  'new_Ca_triple_rule':None,'cost_planning':costs,'new_molecular_calls_in_scope_preparation':0,'production_changed':False}
 out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);write(out/'SCIENTIFIC_SCOPE.json',scope)
 summary={k:scope[k] for k in ('protocol_id','counts','execution_resources','cost_planning')};summary.update(status='prepared_scope_only',scope=pin(out/'SCIENTIFIC_SCOPE.json'),frozen_reference=scope['reference'])
 write(out/'PREFLIGHT.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--output',required=True);a=p.parse_args();prepare(a.config,a.output)
