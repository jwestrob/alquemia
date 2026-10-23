#!/usr/bin/env python3
"""Design, hashes and row-layout metadata only; no scientific recomputation."""
import hashlib,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;P=HERE.parent;PINS={};checks=[]
def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(q):
 assert pin(q['path'])==q,q;PINS[q['path']]=q;return Path(q['path'])
def read(q):return json.loads(verify(q).read_text())
def check(ok,name):
 assert ok,name;checks.append(name)
start,cpu=time.monotonic(),time.process_time()
rp=pin(P/'READY.json');r=read(rp)
check(rp['sha256']=='2c4f7043b556dceba03fc879d68d4e2bffdcf9bdbbfcb2f8bd4b7ae2fcba7bae','Exact frozen READY')
for k in ('root_declaration','design','inventory','inventory_producer','inventory_result'):verify(r[k])
i=read(r['inventory']);receipt=read(r['inventory_result']);dec=read(r['root_declaration'])
check(i['design']==r['design'] and i['declaration']==r['root_declaration'] and i['inventory_implementation']==r['inventory_producer'],'Frozen design/inventory/source joins')
check(receipt['inventory']==r['inventory'] and receipt['implementation']==r['inventory_producer'],'Saved metadata receipt joins')
check(receipt['status']=='metadata_inventory_and_zero_scalar_check_complete' and receipt['negative_inherited_zero_margins']==receipt['negative_new_group_normalized_zero_margins']=={},'Saved zero check has no negatives; not rerun')
check(all(receipt[k]==0 for k in ('new_native_calls','solver_calls','probe_Jacobians_constructed','raw46_replays')),'Saved inventory scope zero native/solver/Jacobian/raw46')
check(r['source_only'] and not r['runnable_solver_implemented'] and r['execution_release'] is None,'Source design only, execution unbound')
check(r['prospective_solves_maximum']==1 and r['target_groups']==6 and r['derivative_vectors']==28,'Single prospective six-group28-probe diagnostic')
check(i['origin_key']==dec['fixed_candidate_key']==r['origin_key']=='5d992dbffb9fbe646102782cb90b53182deff7e384838ffd809683e1a187fa78','Exact final step8 origin')
for k in ('trimer_admission','trimer_result','trimer_inventory','candidate','fit_result','fit_root_admission','proposal','linearization','conic_data','row_layout'):verify(i[k])
for q in i['sources'].values():verify(q)
proposal=read(i['proposal']);layout=read(i['row_layout']);trinv=read(i['trimer_inventory']);admit=read(i['trimer_admission'])
check(proposal['outer_step']==8 and proposal['origin_key']==i['origin_key'] and proposal['linearization']==i['linearization'] and proposal['conic_data']==i['conic_data'] and proposal['row_layout']==i['row_layout'],'Exact saved step8 matrix and provenance')
check(admit['candidate']==i['candidate'] and admit['result']==i['trimer_result'],'Admitted terminal trimer authority')
rows=i['ordered_derivative_vectors']
check(len(rows)==28 and [v['key'] for v in rows]==proposal['derivative_keys'],'Exact ordered28 derivative identity')
for j,v in enumerate(rows):
 verify(v['artifact'])
 check(v['index']==j and v['coordinate']==j//2 and v['side']==('minus' if j%2==0 else 'plus') and not v['raw_native_replay_performed'],'Derivative metadata '+str(j))
 check(set(v['required_dynamic_fields'])==set(i['required_dynamic_cell_ids']) and all(x['status']=='complete' for x in v['required_dynamic_fields'].values()),'Six admitted dynamic fields '+str(j))
check(len(i['actual_spans'])==14 and len(i['required_dynamic_cell_ids'])==6 and len(i['static_terms'])==8,'Finite14-column6+8 input scope')
check(i['cycle_dictionaries']=={k:v for k,v in trinv['cycles'].items() if k.startswith('aquo9_') and k.endswith('__far')},'Literal W9far seven-term cycles')
groups=i['groups'];check([(v['name'],v['dimension']) for v in groups]==[('E3/La',1),('mu3/La',3),('E3/Dy',1),('mu3/Dy',3),('E3/Dy_minus_La',1),('mu3/Dy_minus_La',3)],'All six individual/paired scalar/vector groups')
check(all(v['stored_scale']>0 and v['stored_scale_minus_reconstructed_norm']==0 and v['normalized_zero_margin']==0 for v in groups),'Saved positive scale/zero-margin evidence only')
z=i['zero_step'];check(z['literal_q_ceiling']==i['literal_q_ceiling']==.13804357390849725 and z['q_reconstructed']==.13804357390849722,'Literal cap and saved reconstruction retained')
check(not z['padding_added'] and z['future_actual_conic_zero_check_required'] and z['sqrt_q_squared_minus_literal_q']==0 and z['negative_margins']=={},'No cushion and future constructed-row zero check explicit')
soc=[v for v in layout if v['kind']=='SOC'];lin=[v for v in layout if v['kind']=='nonnegative']
check(len(soc)==20 and soc[0]['start']==0 and soc[-1]['stop']==547 and len(lin)==890 and lin[0]['start']==547 and lin[-1]['stop']==1437,'Original20SOC/890sides/1437rows metadata')
ep=[v for v in layout if v['name'].startswith('epigraph/')]
check([(v['name'],v['start'],v['stop']) for v in ep]==[('epigraph/lower',1433,1434),('epigraph/upper',1434,1435)],'Only two named epigraph rows removed')
check([(v['name'],v['start'],v['stop']) for v in layout[-2:]]==[('signed_W9_z/upper',1435,1436),('signed_W9_z/lower',1436,1437)],'Both signed-W9 sides retained')
check([v['name'] for v in soc]==i['old_SOC_names'],'Original SOC families preserved')
result={'status':'passed_trimer_local_span_design_review','ready':rp,'root_declaration':r['root_declaration'],'design':r['design'],'inventory':r['inventory'],'inventory_producer':r['inventory_producer'],'inventory_result':r['inventory_result'],'reviewer_implementation':pin(__file__),'report':pin(HERE/'REPORT.txt'),'checks':checks,'checks_passed':len(checks),'verified_pins':list(PINS.values()),'verified_pin_count':len(PINS),'expected_matrix_shape':[1455,15],'expected_SOC_count':26,'expected_nonnegative_dimension':890,'exact_arithmetic_equivalence':True,'row_remap':{'old_SOC':[0,547,0,547],'new_SOC':[547,565],'old_linear_before_epigraph':[547,1433,565,1451],'drop_only_old_epigraph':[1433,1435],'old_signed_W9':[1435,1437,1451,1453],'new_t_bounds':[1453,1455]},'zero_step_check_repeated':False,'conic_matrix_or_Jacobian_constructed':False,'executed_solver_calls':0,'new_native_calls':0,'source_mutations':0,'runnable_implementation_admitted':False,'execution_release':False,'blockers':[],'limitations':['This admits the frozen design, not future constructed-row numerical feasibility or solver execution.','Exact q/no-padding and scalar/conic zero checks remain as explicitly declared; failures must be preserved.','Positive feasible t is only a local affine witness; no global model or physical accuracy inference.'],'cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'threads':1}}
with (HERE/'REVIEW.json').open('x') as h:json.dump(result,h,indent=2,allow_nan=False);h.write('\n')
print(json.dumps({'review':pin(HERE/'REVIEW.json'),'report':result['report'],'checks_passed':len(checks),'verified_pin_count':len(PINS),'cost':result['cost']}))
