#!/usr/bin/env python3
"""Independent saved arithmetic; imports no solver or producer implementation."""
import hashlib,json,math,os,sys,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;PINS={};checks=[]
def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(q):
 assert pin(q['path'])==q,q;PINS[q['path']]=q;return Path(q['path'])
def read(q):return json.loads(verify(q).read_text())
def check(ok,name):
 assert ok,name;checks.append(name)
def close(a,b,name):
 x=np.asarray(a,dtype=float);y=np.asarray(b,dtype=float);bound=128*max(24,x.size)*sys.float_info.epsilon*max(1.,float(np.max(abs(x))),float(np.max(abs(y))))
 check(x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all() and float(np.max(abs(x-y)))<=bound,name)
def norm(v):return math.sqrt(math.fsum(float(x)*float(x) for x in v))
def affine(r,J,d):return [float(v)+math.fsum(float(a)*float(b) for a,b in zip(row,d)) for v,row in zip(r,J)]
def load(q):
 with np.load(verify(q),allow_pickle=False) as z:data={k:z[k].copy() for k in z.files}
 a=np.zeros(tuple(data['A_shape']),dtype=float)
 for j in range(a.shape[1]):
  lo,hi=data['A_indptr'][j:j+2];a[data['A_indices'][lo:hi],j]=data['A_data'][lo:hi]
 return a,data['b'],data['c'],data
def cycle_values(v,inv):
 check(v['status']=='complete','Admitted vector status')
 vals={k:(x['energy_kcal_mol'],x['mu_Debye']) for k,x in inv['static_terms'].items()}
 vals.update({k:(v['cells'][k]['energy'],v['cells'][k]['mu_Debye']) for k in inv['required_dynamic_cell_ids']})
 result=[]
 for metal in ('La','Dy'):
  co=inv['cycle_dictionaries']['aquo9_central_'+metal+'__far']['model_terms'];result.extend([[math.fsum(c*vals[k][0] for k,c in co.items())],[math.fsum(c*vals[k][1][j] for k,c in co.items()) for j in range(3)]])
 return result

def main(bind,out):
 tick,cpu=time.monotonic(),time.process_time()
 check(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')),'One thread')
 check(bind['implementation']==pin(__file__) and bind['output_directory']==str(out),'Exact review binding');verify(bind['declaration'])
 r=read(bind['result']);cfg=read(bind['configuration']);proc=read(bind['process_receipt']);ready=read(bind['ready'])
 check(proc['returncode']==0 and proc['result_sha256']==bind['result']['sha256'] and proc['job_id'] is None,'Terminal process receipt')
 check(r['status']=='trimer_local_span_solve_complete' and r['error'] is None and r['mode']=='solve','Completed computation')
 check(r['configuration']==bind['configuration'] and r['implementation']==cfg['implementation']==ready['implementation'],'Executed source/configuration')
 source=read(cfg['source_review']);rel=read(cfg['root_release']);verify(rel['root_source_review'])
 check(source['status']=='passed_trimer_local_span_runtime_source_review' and source['implementation']==cfg['implementation'] and source['configuration_fingerprint']==cfg['configuration_fingerprint']==r['configuration_fingerprint'],'Independent source/fingerprint authority')
 check(rel['status']=='released_trimer_local_span_solve' and rel['mode']=='solve' and rel['maximum_solver_calls']==1,'One exact solve release')
 for k in ('implementation','output_directory','source_review','runtime_declaration','configuration_fingerprint'):check(rel[k]==cfg[k],'Release '+k)
 prep=read(cfg['preparation']);pa=read(rel['preparation_admission'])
 check(pa['status']=='admitted_trimer_local_span_preparation' and pa['result']==cfg['preparation']==rel['preparation'] and prep['implementation']==r['implementation'] and prep['status']=='trimer_local_span_prepared','Exact admitted preparation')
 check(prep['configuration_fingerprint']==r['configuration_fingerprint'] and prep['groups']==pa['groups'],'Same scientific preparation')
 attempt=read(r['attempt']);check(r['attempt']['path']==cfg['attempt_path'] and attempt['configuration_fingerprint']==cfg['configuration_fingerprint'] and attempt['implementation']==r['implementation'] and attempt['maximum_attempts']==1 and attempt['PID']==r['cost']['PID'],'Single durable attempt')
 check(r['solver_calls']==r['solver_constructions']==1 and r['native_calls']==r['new_probes']==0,'One construction/solve; zero new chemistry')
 check(r['cost']['threads']==1 and r['cost']['job_id'] is None and r['cost']['hostname']=='biotite','Actual login resource accounting')
 check(r['solver_status']==r['solution']['status']==r['info']['status'] and r['solver_success']==(r['solver_status']=='Solved'),'Raw solver status unchanged')
 inv=read(cfg['inventory']);p=read(inv['linearization']);groups=read(prep['groups']);oldlayout=read(inv['row_layout']);layout=read(r['row_layout'])
 check(r['inventory']==cfg['inventory']==ready['inventory'] and p['raw_center_q']==inv['literal_q_ceiling'],'Fixed inventory/literal q')
 check(read(inv['trimer_admission'])['result']==inv['trimer_result'],'Prior trimer admission reused')
 settings=read(inv['sources']['settings']);check(r['effective_settings']==prep['effective_settings']==settings,'All effective settings exact')
 check(settings['max_iter']==200 and settings['max_threads']==1 and settings['time_limit']==60. and r['solution']['iterations']<=200,'Frozen finite solver settings')
 verify(cfg['conic_helpers']);verify(inv['sources']['core']);verify(inv['sources']['solver'])
 # Independently sum only the declared narrow fields, not original objective/Jacobian.
 center=read(inv['candidate']);cv=cycle_values(center,inv);probes=[]
 for n,row in enumerate(inv['ordered_derivative_vectors']):
  v=read(row['artifact']);check(v['parameters']==row['parameters'] and row['index']==n,'Exact saved stencil metadata '+str(n));probes.append(cycle_values(v,inv))
 refs=[]
 for metal in ('La','Dy'):
  c=inv['cycle_dictionaries']['aquo9_central_'+metal+'__far'];refs.extend([[c['reference_electronic_E3_kcal_mol']],c['reference_mu3_Debye']])
 gj=[]
 for k in range(4):
  J=[[ (probes[2*j+1][k][a]-probes[2*j][k][a])/span for j,span in enumerate(inv['actual_spans'])] for a in range(len(cv[k]))]
  gj.append({'r':[a-b for a,b in zip(cv[k],refs[k])],'J':J,'model':cv[k],'reference':refs[k]})
 for a,b in ((0,2),(1,3)):
  model=[x-y for x,y in zip(cv[b],cv[a])];ref=[x-y for x,y in zip(refs[b],refs[a])]
  gj.append({'r':[x-y for x,y in zip(model,ref)],'J':(np.asarray(gj[b]['J'])-np.asarray(gj[a]['J'])).tolist(),'model':model,'reference':ref})
 names=['E3/La','mu3/La','E3/Dy','mu3/Dy','E3/Dy_minus_La','mu3/Dy_minus_La']
 check([g['name'] for g in groups]==names,'Six signed scalar/vector groups')
 for g,new,authority in zip(groups,gj,inv['groups']):
  for k in ('r','J','model','reference'):check(np.array_equal(g[k],new[k]),'Exact independent group '+g['name']+'/'+k)
  check(g['r']==authority['r'] and g['scale']==norm(g['r']) and g['scale']>0,'Literal center/scale '+g['name'])
 A,b,c,data=load(r['conic_data']);OA,ob,oc,od=load(inv['conic_data']);PA,pb,pc,pd=load(prep['conic_data'])
 check(all(np.array_equal(data[k],pd[k]) for k in data) and layout==read(prep['row_layout']),'Exact executed/prepared arrays and layout')
 check(A.shape==(1455,15) and OA.shape==(1437,15) and np.array_equal(c,np.r_[np.zeros(14),-1.]),'Final dimensions/objective minus-t')
 oldby={v['name']:v for v in oldlayout};newby={v['name']:v for v in layout};q0=inv['literal_q_ceiling'];s0=math.sqrt(q0)
 check(set(oldby)-set(newby)=={'epigraph/lower','epigraph/upper'},'Only two epigraph rows removed')
 for name,row in newby.items():
  lo,hi=row['start'],row['stop']
  if name in oldby:
   rr=oldby[name];ids=list(range(rr['start'],rr['stop']));check(row['original_rows']==ids and np.array_equal(A[lo:hi,:14],OA[ids,:14]) and np.all(A[lo:hi,14]==0) and np.array_equal(b[lo:hi],ob[ids]-OA[ids,14]*s0),'Exact fixed-s retained row '+name)
  elif name.startswith('trimer/'):
   g=next(g for g in groups if 'trimer/'+g['name']==name);expected=np.zeros((len(g['r'])+1,15));expected[0,14]=1;expected[1:,:14]=-np.asarray(g['J'])/g['scale']
   check(np.array_equal(A[lo:hi],expected) and np.array_equal(b[lo:hi],np.r_[1.,np.asarray(g['r'])/g['scale']]),'Exact new cone '+name)
  else:
   expected=np.zeros(15);expected[-1]=-1 if name=='fraction/lower' else 1
   check(name in ('fraction/lower','fraction/upper') and np.array_equal(A[lo],expected) and b[lo]==(0 if name=='fraction/lower' else 1),'Exact fraction side '+name)
 check(len([v for v in layout if v['kind']=='SOC'])==26 and len([v for v in layout if v['kind']=='nonnegative'])==890,'26SOC/890linear retained denominator')
 # Independent returned scalar arithmetic, never use producer certificate helper.
 x=r['solution']['x'];d=x[:14];t=x[14];m={};blocks={};gates={}
 for name,g in p['blocks'].items():
  v=affine(g['r'],g['J'],d);sq=math.fsum(y*y for y in v);ratio=sq/g['denominator'];blocks[name]={'predicted_squared_norm':sq,'fixed_denominator':g['denominator'],'ratio':ratio}
  if name!='old_EF':m['block/'+name]=(q0 if name in ('P_energy','B_energy') else 1.)-ratio
 for name,v in zip(p['energy']['names'],affine(p['energy']['r'],p['energy']['J'],d)):
  m['energy/upper/'+name]=1-v;m['energy/lower/'+name]=1+v;gates['energy/'+name]={'error_kcal_mol':v,'passed':abs(v)<=1}
 for name,g in p['forces'].items():
  v=affine(g['r'],g['J'],d);rms=math.fsum(float(w)*a*a for w,a in zip(g['weights'],v));maximum=max(abs(a) for a in v);m['force/RMS/'+name]=1-rms/4
  for j,a in enumerate(v):m['force/upper/'+name+'/'+str(j)]=1-a/5;m['force/lower/'+name+'/'+str(j)]=1+a/5
  gates['force/'+name]={'rigid_RMS_kcal_mol_A':math.sqrt(rms),'maximum_body_component_error_kcal_mol_A':maximum,'passed':rms<=4 and maximum<=5}
 signed=p['signed_W9_z'];sv=signed['r_Debye']+math.fsum(float(a)*b for a,b in zip(signed['J_Debye'],d));m.update({'signed_W9/upper':signed['absolute_ceiling_Debye']-sv,'signed_W9/lower':signed['absolute_ceiling_Debye']+sv,'fraction/lower':t,'fraction/upper':1-t})
 params=[v+w*dv for v,w,dv in zip(p['center'],p['width'],d)]
 for j,dv in enumerate(d):
  for name,v in [('trust/lower',dv+.025),('trust/upper',.025-dv),('box/lower',params[j]-p['lower'][j]),('box/upper',p['upper'][j]-params[j])]:m[name+'/'+str(j)]=v
 detail={}
 for g in groups:
  v=affine(g['r'],g['J'],d);nv=norm(v);margin=(1-t)*g['scale']-nv;m['trimer/'+g['name']]=margin;detail[g['name']]={'signed_residual':v,'reference':g['reference'],'predicted_model':[a+b for a,b in zip(v,g['reference'])],'norm':nv,'center_norm':g['scale'],'margin':margin,'achieved_fraction':1-nv/g['scale']}
 q=max(blocks[k]['ratio'] for k in ('P_energy','B_energy'));m['literal_q_cap']=q0-q;feasible=all(v>=0 for v in m.values()) and all(g['passed'] for g in gates.values());wit=r['witness']
 check(m==wit['margins'] and len(m)==917,'All917 returned scalar margins exact')
 check(blocks==wit['blocks'] and gates==wit['original25'] and len(gates)==25 and detail==wit['groups'],'Exact all25/block/six-group scalar results')
 check(wit['parameters']==params and wit['normalized_displacement']==d and wit['fraction']==t and wit['q']==q and wit['signed_W9_residual']==sv,'Exact returned physical point and objective quantities')
 check(wit['minimum_margin']==min(m.values()) and wit['exact_sign_feasible']==feasible and wit['all25_predicted']==all(v['passed'] for v in gates.values()),'No tolerance changes exact feasibility')
 negatives={k:v for k,v in m.items() if v<0}
 # Independent fsum matrix products; tolerances check arithmetic only, never signs.
 raw_s=r['solution']['s'];dual=r['solution']['z'];Ax=[math.fsum(float(a)*float(v) for a,v in zip(row,x)) for row in A];slack=[float(bb)-v for bb,v in zip(b,Ax)]
 primal=[v+ss-float(bb) for v,ss,bb in zip(Ax,raw_s,b)];stationary=[math.fsum(float(A[i,j])*float(dual[i]) for i in range(len(b)))+float(c[j]) for j in range(15)]
 pcost=math.fsum(float(a)*float(v) for a,v in zip(c,x));dcost=-math.fsum(float(a)*float(v) for a,v in zip(b,dual));comp=math.fsum(float(a)*float(v) for a,v in zip(raw_s,dual));diag=r['independent_residuals']
 close(primal,diag['primal_equation_residual'],'Independent raw primal equation residual');close(stationary,diag['dual_stationarity_residual'],'Independent dual stationarity residual')
 for value,key in [(max(map(abs,primal)),'primal_residual_inf'),(max(map(abs,stationary)),'dual_residual_inf'),(pcost,'primal_objective_t'),(dcost,'dual_objective_t'),(pcost-dcost,'signed_gap_t'),(comp,'complementarity_dot')]:close(value,diag[key],'Independent KKT '+key)
 def cone(v):
  ans={}
  for row in layout:
   z=v[row['start']:row['stop']];ans[row['name']]=float(z[0])-norm(z[1:]) if row['kind']=='SOC' else float(z[0])
  return ans
 cone_sets={}
 for name,v in [('primal_affine_cone_margins',slack),('raw_slack_cone_margins',raw_s),('dual_cone_margins',dual)]:
  cm=cone(v);check(set(cm)==set(diag[name]),'Cone diagnostic denominator '+name)
  close(list(cm.values()),list(diag[name].values()),'Independent cone diagnostics '+name);cone_sets[name]={'minimum':min(cm.values()),'negative_margins':{k:v for k,v in cm.items() if v<0}}
 close(pcost,r['solution']['obj_val'],'Raw primal objective is minus fraction');close(dcost,r['solution']['obj_val_dual'],'Raw dual objective arithmetic');close(abs(r['solution']['obj_val']-r['solution']['obj_val_dual']),r['info']['gap_abs'],'Solver gap reports retained')
 check(pcost==-t and r['info']['linsolver']['threads']==1 and r['info']['linsolver']['name']=='qdldl','Objective sign and actual one-thread backend')
 check(r['selected_candidate'] is None and r['native_admitted'] is False and r['global_infeasibility_claim'] is False and wit['native_admitted'] is False,'No native/model/global admission')
 verify(pin(Path(bind['result']['path']).with_name('SOLVER.log')))
 summary={'solver_status':r['solver_status'],'solver_success':r['solver_success'],'fraction_t':t,'minimized_objective_minus_t':pcost,'raw_dual_objective':r['solution']['obj_val_dual'],'raw_gap_abs':r['info']['gap_abs'],'q':q,'literal_q_cap':q0,'q_violation':q-q0,'all25_predicted':wit['all25_predicted'],'exact_sign_feasible':feasible,'negative_scalar_margins':negatives,'group_predictions':detail,'unscaled_primal_residual_inf':max(map(abs,primal)),'unscaled_dual_residual_inf':max(map(abs,stationary)),'reported_scaled_primal_residual':r['info']['res_primal'],'reported_scaled_dual_residual':r['info']['res_dual'],'cone_diagnostics':cone_sets,'interpretation':'Numerical near-boundary affine indication only; no positive feasible witness admitted, no optimum/global family/native claim.' if not feasible else 'Exactly feasible affine witness only; no native/global claim.'}
 result={'status':'passed_trimer_local_span_actual_arithmetic_review','implementation':pin(__file__),'binding':pin(HERE/'CONFIG.json'),'declaration':bind['declaration'],'result':bind['result'],'configuration':bind['configuration'],'process_receipt':bind['process_receipt'],'producer':r['implementation'],'source_review':cfg['source_review'],'root_release':cfg['root_release'],'preparation':cfg['preparation'],'preparation_admission':rel['preparation_admission'],'groups':prep['groups'],'conic_data':r['conic_data'],'row_layout':r['row_layout'],'attempt':r['attempt'],'summary':summary,'scalar_margins':m,'checks':checks,'checks_passed':len(checks),'verified_pins':list(PINS.values()),'verified_pin_count':len(PINS),'source_admission_passed':True,'matrix_mapping_and_narrow_group_arithmetic_passed':True,'returned_scalar_arithmetic_passed':True,'primal_dual_arithmetic_passed':True,'positive_feasible_witness_admitted':feasible and t>0,'solver_convergence_claim':False,'model_admitted':False,'new_native_calls':0,'new_solver_calls':0,'new_reference_calls':0,'original531_Jacobian_replayed':False,'full46_native_review_repeated':False,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'producer_preparation':prep['cost'],'producer_execution':r['cost'],'producer_solve_subset':r['solver_cost'],'inherited_costs_not_added':True}}
 lines=['Saved W9far local-span actual arithmetic review: PASSED','Computation completed one solve; returned status '+r['solver_status']+' remains unchanged.','All917 scalar margins and all25 physical gates reproduced exactly.','t='+repr(t)+'; minimized objective -t='+repr(pcost)+'.','q='+repr(q)+'; exact cap='+repr(q0)+'; violation='+repr(q-q0)+'.','Negative scalar margins: '+json.dumps(negatives,sort_keys=True),'Exact feasible witness: '+str(feasible)+'. No positive witness is admitted when this is false.','Raw primal/dual gap='+repr(r['info']['gap_abs'])+'; unscaled residual infinity norms='+repr([max(map(abs,primal)),max(map(abs,stationary))])+'.','Numerical near-boundary evidence is not certified feasibility, convergence, optimality, global family inadequacy or native accuracy.','No tolerance padding, fallback, solver/native/reference execution or original531/full46 replay.','','Six group achieved fractional reductions:']
 lines += [name+': '+repr(v['achieved_fraction'])+'; margin='+repr(v['margin']) for name,v in detail.items()]
 lines += ['','Actual preparation cost: '+json.dumps(prep['cost'],sort_keys=True),'Actual solve-run cost: '+json.dumps(r['cost'],sort_keys=True),'Solver subset cost (not additional): '+json.dumps(r['solver_cost'],sort_keys=True),'Saved review cost: '+json.dumps(result['cost'],sort_keys=True)]
 (out/'REPORT.txt').write_text('\n'.join(lines)+'\n');result['report']=pin(out/'REPORT.txt');(out/'REVIEW.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'review':pin(out/'REVIEW.json'),'report':result['report'],'summary':{k:v for k,v in summary.items() if k not in ('group_predictions','cone_diagnostics')},'checks_passed':len(checks),'cost':result['cost']}))

if __name__=='__main__':
 tick,cpu=time.monotonic(),time.process_time();bind=read(pin(HERE/'CONFIG.json'));out=Path(bind['output_directory']);out.mkdir(parents=True,exist_ok=False)
 try:main(bind,out)
 except BaseException as exc:
  (out/'FAILED_REVIEW.json').write_text(json.dumps({'status':'failed_trimer_local_span_actual_arithmetic_review','implementation':pin(__file__),'binding':pin(HERE/'CONFIG.json'),'error':{'type':type(exc).__name__,'message':str(exc)},'checks_completed':checks,'new_native_calls':0,'new_solver_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}},indent=2)+'\n');raise
