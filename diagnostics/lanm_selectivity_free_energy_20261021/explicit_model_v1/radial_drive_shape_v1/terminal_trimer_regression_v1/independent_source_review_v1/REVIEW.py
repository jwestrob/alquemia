#!/usr/bin/env python3
"""Saved source/metadata joins only. Never import or execute producer methods."""
import ast
import hashlib
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
P=HERE.parent
verified={}
checks=[]
def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(q):
    assert pin(q['path'])==q,q
    verified[q['path']]=q
    return Path(q['path'])
def read(q):return json.loads(verify(q).read_text())
def check(ok,name):
    assert ok,name
    checks.append(name)
def nested_pins(x):
    if isinstance(x,dict):
        if set(x)=={'path','sha256'}:verify(x)
        else:
            for v in x.values():nested_pins(v)
    elif isinstance(x,list):
        for v in x:nested_pins(v)
def functions(q,names,container=()):
    body=ast.parse(verify(q).read_text()).body
    for name in container:body=next(n for n in body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name).body
    return [next(n for n in body if isinstance(n,ast.FunctionDef) and n.name==name) for name in names]

start,cpu=time.monotonic(),time.process_time()
ready_pin=pin(P/'READY.json');ready=read(ready_pin)
check(ready_pin['sha256']=='c5772fa26d903ca3345a5aa7f2d60cb933ae6b01bdb8ce028f1a392e704fabc7','Exact frozen runtime READY')
nested_pins(ready)
scientific=read(ready['scientific_ready']);inv=read(ready['inventory']);s=inv['sources']
nested_pins(scientific);nested_pins(s)
methods=read(ready['method_inputs']);nested_pins(methods)
decl=read(ready['root_declaration']);parent=read(s['parent_manifest'])
check(read(scientific['source_inputs'])==s and scientific['inventory']==ready['inventory'],'Exact scientific/runtime source inventory')
check(inv['candidate']['artifact']==ready['candidate']==scientific['candidate']==decl['parents']['terminal_candidate'],'Exact terminal candidate authority')
check(inv['candidate']['key']==scientific['candidate_key']==decl['candidate_key']=='5d992dbffb9fbe646102782cb90b53182deff7e384838ffd809683e1a187fa78','Fixed20605 terminal identity')
check(inv['cycles']==parent['cycles'] and len(inv['cycles'])==8,'All eight literal parent cycle dictionaries unchanged')
check(inv['rows']==parent['rows'] and len(inv['rows'])==6,'All six inherited row definitions unchanged')
a,b,c=map(set,(inv['reused_selected'],inv['static_terms'],inv['rows']))
check((len(a),len(b),len(c))==(14,22,6) and not(a&b or a&c or b&c),'Disjoint14selected/22static/6new partition')
check(a|b|c==set().union(*(set(v['model_terms']) for v in inv['cycles'].values())),'Complete42-term dependency set')
check(all(f'aquo9_central_{m}__far__metal_two_water' in a for m in ('La','Dy')),'Both W9far selected trimers reused')
check(len(inv['contrasts'])==4 and all(v['sign']=='Dy-minus-La' for v in inv['contrasts'].values()),'Four signed Dy-minus-La contrasts')
check(all(not r['needs_export'] and not r['needs_gradient'] for r in inv['rows'].values()),'Fresh rows cannot reach export/gradient branches')
check(inv['historical_six_row_import_authorized'] is False,'Historical12 six-row import excluded')
check(ready['counts']==inv['counts']==scientific['counts'] and ready['counts']['maximum_new_analyze']==6,'Exact declared execution denominators')
check(ready['native_induction_blocks_expected']==12 and ready['new_optimizer_or_conic_solves']==0,'Internal induction distinct from optimizer calls')
check(inv['physical_accuracy_qualified'] is False and inv['execution_release'] is None,'Scientific accuracy and runtime remain unadmitted')
check(inv['reference']['historical_reference']==parent['historical_reference'],'Original reference and precision metadata copied unchanged')
inventory_result=read(scientific['metadata_result'])
check(inventory_result['status']=='prepared_radial14_terminal_trimer_inventory' and len(inventory_result['checks'])==53 and inventory_result['new_native_calls']==0 and inventory_result['scalar_replays']==0,'Saved53 metadata checks without rerun')
controls=read(ready['controls'])
check(controls['status']=='passed_radial14_terminal_trimer_source_controls' and controls['runner']==ready['implementation'] and controls['implementation']==ready['controls_source'],'Final source control identity')
check(controls['checks_passed']==len(controls['checks'])==15 and all(x['passed'] for x in controls['checks']),'Saved15 manufactured controls passed')
check(controls['actual_native_calls']==controls['actual_optimizer_calls']==controls['actual_scientific_replays']==0,'Controls used no scientific executions')
pre=read(ready['metadata_preflight'])
check(pre['status']=='passed_radial14_terminal_trimer_metadata_preflight' and pre['implementation']==ready['implementation'] and pre['error'] is None,'Final-source preflight passed')
check(pre['native_program_invocations']==0 and pre['cycles']=={} and pre['comparison']=={},'Preflight no native or cycle calculation')
check(pre['configuration']['_configuration_file']==ready['config_template'],'Preflight exact config template')
cfg=read(ready['config_template']);release=read(ready['release_template'])
check(all(cfg[k]==ready[k] for k in ('implementation','scientific_ready','inventory','method_inputs','candidate')),'Configuration joins exact frozen scientific and execution sources')
check(cfg['source_review'] is None and cfg['root_release'] is None and cfg['output_directory'] is None,'Actual config deliberately unbound')
check(release['status']=='UNRELEASED' and release['maximum_native_calls']==6 and release['native_workers']==1 and release['maximum_aggregate_numerical_threads']==2 and release['aggregate_login_limit']==16,'Unreleased finite serial-login resource contract')
check(release['active_process_check'] is None and release['execution_mode']=='login','Future actual capacity receipt required')
runner=verify(ready['implementation']).read_text();tree=ast.parse(runner)
check("extract(s['parent_adapter'],('cycles','comparisons'),{'require':require,'math':math})" in runner,'Exact inherited scalar function extraction')
check("extract(s['radial_adapter'],('checked_values','render','_parse'),radial_ns,('evaluator_class','RadialEvaluator'))" in runner,'Exact qualified radial method extraction')
check("extract(extra['radial_checks'],('gate','delta','readback'),ns)" in runner,'Exact qualified radial readback extraction')
functions(s['parent_adapter'],('cycles','comparisons'))
functions(s['radial_adapter'],('checked_values','render','_parse'),('evaluator_class','RadialEvaluator'))
functions(methods['radial_checks'],('gate','delta','readback'))
check(not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in ('DefaultSolver','solve','evaluate','propose') for n in ast.walk(tree)),'Runner has no optimizer/proposal/full-vector call')
check("old_component_terminal':'feca12_physical_gate_terminal'" in runner and "new_terminal':'frozen_radial14_terminal5d992'" in runner,'Explicit predecessor and candidate labels')
check("result['native_program_invocations']+=1;persist()" in runner and "result['native_program_invocations']<6" in runner,'Durable bounded invocation intent')
check("not os.environ.get('SLURM_JOB_ID')" in runner and "not os.environ.get('SLURM_STEP_ID')" in runner and "'controller_pid':os.getpid()" in runner,'Truthful login host/PID/no fake scheduler')
wrapper=verify(ready['wrapper']).read_text()
check(all(k+'=1' in wrapper for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')) and 'SLURM' not in wrapper,'Wrapper sets six one-thread limits without erasing scheduler context')
check(ready['native_release'] is None and ready['new_native_calls']==0,'Frozen source is not a runtime release')
result=dict(ready['required_source_review'])
result.update(ready=ready_pin,scientific_specification=scientific['scientific_specification'],root_declaration=ready['root_declaration'],candidate=ready['candidate'],controls=ready['controls'],metadata_preflight=ready['metadata_preflight'],inventory_metadata_result=scientific['metadata_result'],reviewer_implementation=pin(__file__),report=pin(HERE/'REPORT.txt'),checks=checks,checks_passed=len(checks),verified_pins=list(verified.values()),verified_pin_count=len(verified),blockers=[],scope='Independent source and saved metadata review only; no author controls/preparation/scientific execution repeated',new_native_calls=0,new_optimizer_calls=0,new_scientific_replays=0,source_mutations=0,actual_candidate_admission=False,execution_release=False,cost={'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'threads':1})
with (HERE/'REVIEW.json').open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
print(json.dumps({'review':pin(HERE/'REVIEW.json'),'checks_passed':len(checks),'verified_pin_count':len(verified),'cost':result['cost']}))
