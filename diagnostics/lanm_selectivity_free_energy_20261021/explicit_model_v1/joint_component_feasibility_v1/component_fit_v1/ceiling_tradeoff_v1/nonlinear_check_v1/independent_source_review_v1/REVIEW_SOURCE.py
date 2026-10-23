"""One static pin/AST review; no producer imports or author-control execution."""
import ast
import hashlib
import json
from pathlib import Path
import time

tick,cpu=time.monotonic(),time.process_time()
here=Path(__file__).resolve().parent;source_dir=here.parent
root=next(p for p in here.parents if (p/'diagnostics').is_dir() and (p/'workspaces').is_dir())
relative=here.relative_to(root/'diagnostics')
out=root/'workspaces'/relative;out.mkdir(parents=True,exist_ok=False)
checks=[];pins=[]
def require(ok,label):
    if not ok:raise ValueError(label)
    checks.append(label)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(ref):
    require(pin(ref['path'])==ref,'Pin '+ref['path']);pins.append(ref)
    return json.loads(Path(ref['path']).read_text())
def verify(ref):
    require(pin(ref['path'])==ref,'Pin '+ref['path']);pins.append(ref)
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

ready=read(pin(source_dir/'READY.json'))
require(pin(source_dir/'READY.json')['sha256']=='bdf3863e93d432408b5dd4a68f5254da985555470de62e525d55c51f1796d03d','Exact assigned frozen readiness')
for v in ready.values():
    if isinstance(v,dict) and set(v)=={'path','sha256'}:verify(v)
cfg=read(ready['configuration']);template=read(ready['run_template']);prep=read(ready['preparation']);controls=read(ready['controls'])
fingerprint=lambda c:digest({k:v for k,v in c.items() if k not in ('root_release','source_review','preparation','tradeoff_actual_review')})
require(fingerprint(cfg)==fingerprint(template)==ready['configuration_fingerprint']==prep['configuration_fingerprint'],'Same scientific preparation/run fingerprint')
for k in ('parent_evaluator','parent_configuration','followup_policy','tradeoff_result','tradeoff_actual_review','original_pilot_implementation'):verify(cfg[k])
review=read(cfg['tradeoff_actual_review']);declaration=read(cfg['declaration']);policy=read(cfg['followup_policy'])
candidate=cfg['candidate']
require(candidate==ready['candidate']==template['candidate']==prep['candidate']==declaration['candidate'],'Frozen candidate in every binding')
require(candidate['fraction']==1. and candidate['key']==digest(candidate['parameters'])=='d1fc95754c0cd8dceb9a7616caf2e8409974f6952d0a904e0f2ead9645acbec3','Exact full physical candidate')
require(review['status']=='passed_ceiling_tradeoff_actual_arithmetic_review' and review['certificate_checks_passed'] is True and review['result']==cfg['tradeoff_result'] and review['prospective_policy']==cfg['followup_policy'],'Actual saved arithmetic/policy joins')
require(review['prospective_trigger_witness']=={'fraction':1.,'physical_key':candidate['key'],'parameters':candidate['parameters']} and review['solver_status']=='AlmostSolved' and review['solver_success'] is False,'Unchanged reviewed witness and reduced solver status')
require(policy['selection_rule']['order']==[1.,.5,.25,.125] and policy['selection_rule']['minimum_absolute_gain_in_component_q']==.05,'Prospective order and trigger')
old=ast.parse(Path(cfg['original_pilot_implementation']['path']).read_text());new=ast.parse(Path(ready['implementation']['path']).read_text())
def functions(tree):return {n.name:n for n in ast.walk(tree) if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
of,nf=functions(old),functions(new)
for name in ('allowed','_stop','_invoke','_cell','evaluate','snapshot','one','__init__'):
    require(ast.dump(of[name],include_attributes=False)==ast.dump(nf[name],include_attributes=False),'Identical inherited pilot AST '+name)
require(all(n not in nf for n in ('render','_parse','_assemble','cache_identity')),'Scientific evaluator methods inherited without override')
parent=read(cfg['parent_configuration']);require(parent['implementation']==cfg['parent_evaluator'],'Exact parent evaluator')
manifest=read(parent['manifest']);require(set(prep['prepared_cells'])==set(manifest['rows']) and len(manifest['rows'])==46,'Prepared exact46 manifest cells')
require(sum(v['needs_gradient'] for v in manifest['rows'].values())==10 and sum(v['needs_export'] for v in manifest['rows'].values())==4,'46analyze+10testgrad/four export fixtures')
require(prep['status']=='prepared_ceiling_tradeoff_native_zero_calls' and prep['error'] is None and prep['native_programs']==0 and prep['snapshot']['native_program_invocations']==0,'Actual zero-call preparation only')
require(controls['status']=='passed_ceiling_tradeoff_native_source_controls' and controls['error'] is None and controls['check_count']==len(controls['checks'])==23 and controls['implementation']==ready['implementation'],'Existing23 source controls inspected, not rerun')
require(cfg['limits']=={'new_parameter_vectors':1,'native_program_attempts':56,'workers':12,'aggregate_login_thread_budget':16},'Finite execution caps')
require(cfg['executor']=={'mode':'login','hostname':'biotite','workers':12,'numerical_threads':1,'aggregate_thread_budget':16},'12native+1controller=13 within16 aggregate')
qualification=read(parent['qualification_review']);build=read(parent['build_receipt'])
require(qualification['status']=='passed_local_response_native_qualification' and all(qualification[k] is True for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed')),'Existing native numerical admission')
require(all(qualification[k]==parent[k] for k in ('source_manifest','build_receipt','executables')) and build['manifest']==parent['source_manifest'],'Exact source/build/binary lineage')
verify(parent['objective'])
report={'status':'passed_ceiling_tradeoff_native_source_review','implementation':ready['implementation'],
 'configuration_fingerprint':ready['configuration_fingerprint'],'reviewer_source':pin(__file__),
 'readiness':pin(source_dir/'READY.json'),'configuration':ready['configuration'],'preparation':ready['preparation'],
 'controls':ready['controls'],'actual_tradeoff_review':cfg['tradeoff_actual_review'],'followup_policy':cfg['followup_policy'],
 'candidate':candidate,'parent_configuration':cfg['parent_configuration'],
 'authorship_disclosure':'Reviewer performed the saved tradeoff arithmetic review and authored the earlier conic adapter. This nonlinear pilot source was authored by Arendt. No producer or test function was executed here.',
 'findings':[
  'Only the reviewed full-fraction physical vector is permitted. Auxiliary q never enters native parameters. Prospective policy and actual review are bound at the first native invocation, preserving AlmostSolved.',
  'Invocation, per-cell execution, allowlist, initialization, stop/accounting, evaluate, snapshot and literal C/P/B-cycle bodies match the admitted prior pilot AST. Rendering, parser, cache identity and assembly are inherited unchanged from the exact qualified evaluator.',
  'Fresh immutable output/cache; one parameter-vector request;56 invocation-attempt cap;46 complete analyze cells plus10 testgrad calls required. Global first-failure stop prevents further invocation boundaries; already-started work is accounted, failed/partial/untouched rows retained. No retry or fallback fraction.',
  '12 single-thread native workers plus one numerical controller conservatively use13 numerical threads, within the16 aggregate login budget. Thread environment, hostname and null scheduler identity are checked; root remains responsible for aggregate concurrent ownership.',
  'All25 original energy/force gates remain computed without changed thresholds. Completion is a complete experimental report, not gate/model admission. Old four score-ceiling losses remain descriptive and do not reject a completed batch.',
  'Literal paired C/P/B errors retain kcal/mol and J_P/J_B retain squared kcal/mol. All original block scores and changes are saved. Square roots of the objective-defined pair/C/chi blocks produce equal-parent Cartesian RMS Debye, four-parent Cartesian RMS Debye and four-state RMS a0^3 respectively. Predicted old block norms are multiplied by4 to undo the original531 half-residual weighting before the square root.',
  'Runtime review/preparation/release bindings may be filled without changing the scientific fingerprint. Actual tradeoff review is excluded only from preparation fingerprint and must match exact result/policy/status/witness in the runtime gate and root release.'
 ],'checks':checks,'pin_occurrences':len(pins),'new_native_calls':0,'new_solver_calls':0,'author_controls_rerun':False,'preparation_rerun':False,
 'limits':'Source/preparation readiness only. Native outcome, nonlinear transfer, response accuracy and model/MD admission remain unavailable.',
 'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}}
dest=out/'REVIEW.json';dest.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'review':pin(dest),'checks':len(checks),'pin_occurrences':len(pins)}))
