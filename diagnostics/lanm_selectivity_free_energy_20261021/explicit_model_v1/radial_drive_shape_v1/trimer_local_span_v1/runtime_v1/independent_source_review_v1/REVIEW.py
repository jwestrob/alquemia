#!/usr/bin/env python3
"""Source/pin and saved control joins only; never import runtime numerical code."""
import ast,hashlib,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;P=HERE.parent;PINS={};checks=[]
def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(q):
 assert pin(q['path'])==q,q;PINS[q['path']]=q;return Path(q['path'])
def read(q):return json.loads(verify(q).read_text())
def check(ok,name):
 assert ok,name;checks.append(name)
def nested(x):
 if isinstance(x,dict):
  if set(x)=={'path','sha256'}:verify(x)
  else:
   for v in x.values():nested(v)
 elif isinstance(x,list):
  for v in x:nested(v)
start,cpu=time.monotonic(),time.process_time();rp=pin(P/'READY.json');r=read(rp)
check(rp['sha256']=='e0ef6ce2cc6b591665f3e90d6ee3162d27cad28311bb342ea3abc8e61d50e107','Exact frozen runtime READY');nested(r)
cfg=read(r['configuration_template']);inv=read(r['inventory']);controls=read(r['controls']);decl=read(r['runtime_declaration']);design=read(r['design_admission'])
for q in inv['sources'].values():verify(q)
nested(cfg)
fpkeys=('implementation','inventory','design_ready','conic_helpers','attempt_path','limits','runtime_declaration','design_admission')
fp=hashlib.sha256(json.dumps({k:cfg[k] for k in fpkeys},sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
check(fp==cfg['configuration_fingerprint']==r['configuration_fingerprint']=='2dd73987c1c2b8400f92ed0689d1ce1fa351dcff30df7422ea48ca37c72230e0','Exact acyclic scientific fingerprint')
check(all(cfg[k]==r[k] for k in ('implementation','inventory','design_ready','runtime_declaration','design_admission')),'Source/input declaration joins')
check(decl['input_authority']==cfg['inventory'] and decl['design_admission']==cfg['design_admission'] and design['inventory']==cfg['inventory'] and design['ready']==cfg['design_ready'],'Exact root design/runtime authority')
check(cfg['limits']=={'solver_calls':1,'native_calls':0,'new_probes':0,'threads':1,'solver_time_limit_seconds':60},'Finite one-solve scope')
check(all(cfg[k] is None for k in ('source_review','root_release','preparation','output_directory')),'Templates remain unreleased')
check(controls['status']=='passed_manufactured_runtime_controls' and controls['implementation']==r['implementation'] and controls['test_source']==r['control_source'],'Saved controls bind final source')
check(controls['count']==len(controls['checks'])==28,'All28 saved manufactured controls passed')
check(all(controls[k]==0 for k in ('actual_solver_calls','actual_stencil_reads','actual_Jacobians_constructed','native_calls')),'No actual scientific execution in controls')
source=verify(r['implementation']).read_text();tree=ast.parse(source)
defs={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
check(set(('build_groups','adapt','certificate','gate','release_gate','once','run'))<=set(defs),'All reviewed bounded entrypoints present')
check("edits=ns['clone'](module,ns,'certify_fraction',{13:15,12:14})" in source,'Only scalar checker dimension edits')
check("extract(inv['sources']['original_cycle_helper'],['cycles'],ns)" in source,'Original cycle arithmetic reused')
check("(np.asarray(probes[2*j+1][i])-np.asarray(probes[2*j][i]))/span" in source and "inv['actual_spans']" in source,'Saved plus-minus/span secants')
check("np.array_equal(np.asarray(g['r']),np.asarray(authority['r']))" in source,'Literal center residual equality before scaling')
check("v=b[ids]-np.asarray(A[ids,14].toarray()).ravel()*fixed" in source,'Correct fixed-s RHS substitution')
check("a[0,14]=1.;a[1:,:14]=-J/s" in source,'Correct six-cone slack signs')
check("removed==['epigraph/lower','epigraph/upper']" in source,'Only two original rows dropped')
check("lower[0,14]=-1" in source and "add('fraction/upper','nonnegative',-lower,[1.])" in source,'Exact t lower/upper signs')
check("c=np.r_[np.zeros(14),-1.]" in source,'Maximize t via objective minus-t')
check("old=checker(p,d+[q0],1.,False)['certificate']" in source and "if k!='block/old_EF'" in source,'Literal q scalar certificate; only obsolete block excluded')
check("m['literal_q_cap']=q0-q" in source and "exact_sign_feasible=all(v>=0 for v in m.values()) and old['all25_predicted']" in source,'Exact nonrelaxed scalar signs and original25')
check("require(zero['exact_sign_feasible'] and min(zm.values())>=0" in source,'Constructed cone and scalar zero checks precede solve')
check("np.array_equal(expected[k],v)" in source and "read(prep['row_layout'])==model['layout']" in source,'Exact successful preparation arrays/layout reused')
check("with lock.open('x')" in source and 'os.fsync(f.fileno())' in source,'Exclusive durable attempt before construction')
check("status='attempt_reserved',attempt=pin(lock),solver_constructions=1" in source and "status='solve_entered',solver_calls=1" in source,'Constructor/solve counters retained separately')
check(sum(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='solve' for n in ast.walk(tree))==1,'Only one explicit solve site')
check(not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in ('invoke','evaluate','propose','Popen') for n in ast.walk(tree)),'No native/vector/fit invocation site')
check("solver_status=str(sol.status),solver_success=str(sol.status)=='Solved'" in source and 'independent_residuals=' in source,'Raw status and numerical residuals retained')
check("status='failed_trimer_local_span_no_retry'" in source and 'finally:' in source,'Failed terminal status and final cost persistence')
check("selected_candidate=None,native_admitted=False,global_infeasibility_claim=False" in source,'No candidate or global model admission')
commands=read(r['commands']);check(all(v['environment'][k]=='1' for v in commands.values() for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')),'One-thread commands for both modes')
check(not r['actual_saved_stencil_preparation_performed'] and r['actual_solver_calls']==r['native_calls']==0 and r['execution_release'] is None,'Source readiness distinct from actual release')
result=dict(r['source_review_contract']);result.update(ready=rp,specification=r['specification'],configuration_template=r['configuration_template'],controls=r['controls'],control_source=r['control_source'],reviewer_implementation=pin(__file__),report=pin(HERE/'REPORT.txt'),checks=checks,checks_passed=len(checks),verified_pins=list(PINS.values()),verified_pin_count=len(PINS),blockers=[],source_edits_requested=False,saved_controls_rerun=False,actual_stencil_preparation_performed=False,new_Jacobians_constructed=0,new_solver_calls=0,new_native_calls=0,execution_release=False,scope='Independent source and saved control metadata review only',reporting_note='Inherited primal_objective_t/dual_objective_t report minimized -t; actual common improvement is solution.x[14]/witness.fraction.',cost={'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'threads':1})
with (HERE/'REVIEW.json').open('x') as h:json.dump(result,h,indent=2,allow_nan=False);h.write('\n')
print(json.dumps({'review':pin(HERE/'REVIEW.json'),'report':result['report'],'checks_passed':len(checks),'verified_pin_count':len(PINS),'cost':result['cost']}))
