#!/usr/bin/env python3
"""Bind the unchanged admitted 57-input preparation to the new vacuum binary.

Metadata/source checks only. No baseline numerical re-admission or native call.
"""
from pathlib import Path
import ast
import copy
import difflib
import hashlib
import json
import time

HERE = Path(__file__).resolve().parent
PROTO = HERE.parent
ROOT = PROTO.parents[5]
J = 'explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1'
D = ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/J
W = ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/J
PARENT = D/'sampling_cpu_refactor_v1/qualification'
OUT = W/'periodic_native_prototype_v1/qualification_preparation_v1'


def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def main():
    wall=time.monotonic();cpu=time.process_time();checks=[]
    def check(name, ok):
        checks.append({'name':name,'passed':bool(ok)})
        if not ok: raise AssertionError(name)
    parent_ready_path=D/'sampling_cpu_refactor_v1/LOGIN_QUALIFICATION_READY.json'
    ready=json.loads(parent_ready_path.read_text())
    parent_preparation=json.loads(Path(ready['manifest']['path']).read_text())
    check('parent_manifest_pin',pin(ready['manifest']['path'])==ready['manifest'])
    check('parent_source_pin',pin(ready['implementation']['path'])==ready['implementation'])
    admission_path=W/'sampling_cpu_refactor_v1/qualification_login_run_v1/ROOT_REVIEW.json'
    admission=json.loads(admission_path.read_text())
    check('parent_actual_vacuum_admission',
          admission['status']=='admitted_cpu_refactor_vacuum_numerical_equivalence')
    parent_source=Path(ready['implementation']['path']).read_text()
    old="audit['status']=='passed_local_response_native_source_review'"
    new="audit['status']=='passed_periodic_native_prototype_source_review'"
    check('single_source_audit_status_join',parent_source.count(old)==1)
    source=parent_source.replace(old,new)
    p=HERE/'RUN_LOGIN.py';check('new_runner_not_preexisting',not p.exists());p.write_text(source)
    (HERE/'SOURCE_REVIEW_JOIN.patch').write_text(''.join(difflib.unified_diff(
        parent_source.splitlines(True),source.splitlines(True),
        fromfile=ready['implementation']['path'],tofile=str(p))))
    parent_functions={n.name:ast.dump(n,include_attributes=False)
                      for n in ast.parse(parent_source).body if isinstance(n,ast.FunctionDef)}
    functions={n.name:ast.dump(n,include_attributes=False)
               for n in ast.parse(source).body if isinstance(n,ast.FunctionDef)}
    check('all_scientific_and_execution_functions_AST_identical',all(
          functions[k]==v for k,v in parent_functions.items() if k!='context'))
    check('context_only_declared_status_name_change',source.replace(new,old)==parent_source)
    plan_prefix='''# Vacuum regression for the isolated periodic source prototype

This reuses exactly the admitted parent 57-program/636-internal-FD panel below,
with the same nine systems, input bytes, profiles, methods, parsers and bounds.
Every native program must run anew against the new analyze/testgrad binaries;
old outputs are baseline authorities only, never replacements for new slots.
No periodic/BIO program or candidate composition is included. Both periodic
fatal hooks remain present. Compilation alone does not admit this panel.

The execution envelope here is serial login (one native worker; one numerical
thread per process), replacing only the inherited plan's prospective 12-CPU
Slurm paragraph. Root must review the completed build and this source, then
bind a fresh configuration with root_authorized=true before any execution.
The parent staged first-failure stop and every original check remain intact.

## Unchanged inherited scientific panel

'''
    (HERE/'QUALIFICATION_PLAN.md').write_text(plan_prefix+(PARENT/'QUALIFICATION_PLAN.md').read_text())
    (HERE/'SOURCE_INPUTS.json').write_bytes((PARENT/'SOURCE_INPUTS.json').read_bytes())
    for k,v in parent_preparation['sources'].items():
        check('helper_authority_'+k,pin(v['path'])==v)
    programs=parent_preparation['programs'];rows=parent_preparation['rows']
    check('exact_57_unique_slots',len(programs)==57 and len({p['id'] for p in programs})==57)
    check('exact_636_internal_FD',sum(6*rows[p['cell']]['atoms'] for p in programs if p['fd'])==636)
    check('nine_fixtures',len(rows)==9)
    for job in programs:
        check('exact_three_inputs_'+job['id'],all(pin(v['path'])==v for v in job['input'].values()))
    OUT.mkdir(parents=True,exist_ok=False)
    m=copy.deepcopy(parent_preparation)
    m.update(implementation=pin(p),plan=pin(HERE/'QUALIFICATION_PLAN.md'),
             inherited_preparation=ready['manifest'],inherited_source=ready['implementation'],
             inherited_actual_admission=pin(admission_path),
             preparation_implementation=pin(__file__),
             cost={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,
                   'numerical_threads':1,'new_native_calls':0})
    check('all_program_and_row_records_exactly_inherited',
          m['programs']==parent_preparation['programs'] and m['rows']==parent_preparation['rows'])
    (OUT/'MANIFEST.json').write_text(json.dumps(m,indent=2)+'\n')
    manifest_pin=pin(OUT/'MANIFEST.json')
    source_pin=pin(W/'periodic_native_prototype_v1/source_v1/MANIFEST.json')
    envelope={'status':'declared_serial_vacuum_regression_unreleased',
      'implementation':pin(p),'source_manifest':source_pin,
      'execution':{'mode':'login','hostname':'biotite','job_id':None,'native_workers':1,
        'numerical_threads_per_process':1,'controller_numerical_threads':1,
        'aggregate_login_limit':16,'root_live_resource_review_required':True},
      'native_programs':57,'internal_FD_energies':636,'old_native_slot_reuse':0,
      'periodic_execution_enabled':False,'execution_released':False}
    (HERE/'QUALIFICATION_LOGIN_CONTEXT.json').write_text(json.dumps(envelope,indent=2)+'\n')
    cfg=copy.deepcopy(json.loads((W/'sampling_cpu_refactor_v1/qualification_login_preparation_v1/CONFIG_RUN_UNBOUND.json').read_text()))
    cfg.update(implementation=pin(p),manifest=manifest_pin,plan=m['plan'],
        analyze=None,testgrad=None,source_manifest=source_pin,build_receipt=None,
        source_audit=pin(PROTO/'ROOT_SOURCE_REVIEW.json'),root_authorized=False,
        resource_context=pin(HERE/'QUALIFICATION_LOGIN_CONTEXT.json'))
    (OUT/'CONFIG_RUN_UNBOUND.json').write_text(json.dumps(cfg,indent=2)+'\n')
    result={'status':'prepared_exact_inherited_vacuum57_no_native_calls',
      'implementation':pin(__file__),'runner':pin(p),'manifest':manifest_pin,
      'checks':checks,'check_count':len(checks),'inherited_ready':pin(parent_ready_path),
      'inherited_preparation':ready['manifest'],'new_native_calls':0,'new_numerical_replay':False,
      'cost':{'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    (OUT/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    ready_out={'status':'prepared_vacuum57_build_and_release_binding_pending',
      'implementation':pin(p),'manifest':manifest_pin,'preparation':pin(OUT/'RESULT.json'),
      'configuration_template':pin(OUT/'CONFIG_RUN_UNBOUND.json'),'plan':m['plan'],
      'resource_context':pin(HERE/'QUALIFICATION_LOGIN_CONTEXT.json'),
      'source_manifest':source_pin,'source_audit':cfg['source_audit'],
      'native_programs':57,'internal_FD_energies':636,'native_workers':1,
      'native_slot_reuse':0,'new_native_calls':0,'execution_released':False,
      'bind_after_build':['build_receipt','analyze','testgrad'],
      'bind_after_root_review':['fresh run configuration','root_authorized=true','root release'],
      'argv':['/home/jwestrob/.pyenv/versions/3.12.0/bin/python3',str(p),'run',
              '--config','<ROOT_BOUND_CONFIG>','--out','<FRESH_EXECUTION>'],
      'environment':{k:'1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS',
                  'MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','PYTHONDONTWRITEBYTECODE')}}
    (HERE/'READY.json').write_text(json.dumps(ready_out,indent=2)+'\n')
    print(json.dumps({'ready':pin(HERE/'READY.json'),'result':pin(OUT/'RESULT.json'),
                      'checks':len(checks),'new_native_calls':0}))


if __name__=='__main__':main()
