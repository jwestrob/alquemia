#!/usr/bin/env python3
"""Prepare the existing saved-only actual/parity reviewers; never execute them."""
from pathlib import Path
import ast
import difflib
import hashlib
import json
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[5]
J=Path('explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1')
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/J
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/J
PARENT=D/'sampling_cpu_refactor_v1'
WORK=W/'periodic_native_prototype_v1'


def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def replace(s, old, new):
    assert s.count(old)==1,(old,s.count(old))
    return s.replace(old,new)


def main():
    wall=time.monotonic();cpu=time.process_time()
    result_path=WORK/'qualification_run_v1/execution/RESULT.json'
    # No partial-result interpretation; complete producer is a required bind.
    result=json.loads(result_path.read_text())
    assert result['status']=='local_response_numerical_qualification_complete'
    assert result['error'] is None and result['numerical_checks_passed'] is True
    cfg_path=WORK/'qualification_run_v1/CONFIG_RUN.json'
    cfg=json.loads(cfg_path.read_text());assert cfg==result['config']
    release_path=WORK/'qualification_run_v1/ROOT_RELEASE.json'
    release=json.loads(release_path.read_text())
    assert release['status']=='released_exact_vacuum57_numerical_regression'
    assert release['configuration']==pin(cfg_path)
    ad=HERE/'independent_actual_review_v1';pd=HERE/'independent_parity_review_v1'
    ad.mkdir(exist_ok=False);pd.mkdir(exist_ok=False)
    oldfile=PARENT/'independent_actual_review_v1/REVIEW_ACTUAL.py'
    old=oldfile.read_text();source=old
    source=replace(source,"joint_repair_design_v1/sampling_cpu_refactor_v1')",
                   "joint_repair_design_v1/periodic_native_prototype_v1')")
    oldjoin="    login=q.read(binding['login_review']);ck(login['status']=='passed_sampling_cpu_serial_qualification_source_review' and login['implementation']==binding['runner'] and login['source_manifest']==binding['source_manifest'] and login['resource_context']==cfg['resource_context'],'independent serial-login source review')"
    newjoin="    release=q.read(binding['execution_release']);ck(release['status']=='released_exact_vacuum57_numerical_regression' and release['configuration']==binding['configuration'] and release['source_manifest']==binding['source_manifest'] and release['source_review']==binding['source_review'],'exact reviewed serial-vacuum execution release')"
    source=replace(source,oldjoin,newjoin)
    source=replace(source,"login_source_review=binding['login_review']",
                   "execution_release=binding['execution_release']")
    (ad/'REVIEW_ACTUAL.py').write_text(source)
    (ad/'REVIEWER_BINDINGS.patch').write_text(''.join(difflib.unified_diff(
        old.splitlines(True),source.splitlines(True),fromfile=str(oldfile),
        tofile=str(ad/'REVIEW_ACTUAL.py'))))
    # Reversibility proves every numerical/raw comparison is unchanged.
    restored=source.replace(newjoin,oldjoin).replace(
        "execution_release=binding['execution_release']",
        "login_source_review=binding['login_review']").replace(
        "joint_repair_design_v1/periodic_native_prototype_v1')",
        "joint_repair_design_v1/sampling_cpu_refactor_v1')")
    assert restored==old
    ast.parse(source)
    spec=(PARENT/'independent_actual_review_v1/SPEC.md').read_text()
    spec=spec.replace('build/executable and both source-review joins',
                      'build/executable, root source-review and exact execution-release joins')
    spec += '''
## This source adaptation

The existing independent numerical/raw-output algorithm is unchanged. Only the
prototype path and provenance join change: the actual root source review and
the exact reviewed 57-program execution release replace the parent's separate
serial-resource review. That release is authorization provenance, not an
additional independent scientific review. The source-prototype author prepared
this adapter; root must assign its independent execution/admission. This
preparation performs no raw numerical replay and claims no actual review pass.
'''
    (ad/'SPEC.md').write_text(spec)
    binding={'runner':cfg['implementation'],'source_manifest':cfg['source_manifest'],
             'source_review':cfg['source_audit'],'execution_release':pin(release_path),
             'result':pin(result_path),'configuration':pin(cfg_path),
             'build_receipt':cfg['build_receipt']}
    for p in binding.values():assert pin(p['path'])==p
    (ad/'CONFIG_RUN.json').write_text(json.dumps(binding,indent=2)+'\n')
    authors={'parent_reviewer':pin(oldfile),
             'parent_review':pin(W/'sampling_cpu_refactor_v1/independent_actual_review_v1/actual_v1/REVIEW.json'),
             'preparer':pin(__file__),'methodology_unchanged':True,
             'prepared_by_prototype_author':True,'review_execution_delegated':True}
    (ad/'SOURCE_INPUTS.json').write_text(json.dumps(authors,indent=2)+'\n')
    compare=PARENT/'independent_parity_review_v1/COMPARE.py'
    (pd/'COMPARE.py').write_bytes(compare.read_bytes())
    parity_inputs=json.loads((PARENT/'independent_parity_review_v1/INPUTS.json').read_text())
    parity_inputs.update(new_build=cfg['build_receipt'],new_implementation=cfg['implementation'],
                        new_source=cfg['source_manifest'],expected_new_result=str(result_path))
    (pd/'INPUTS.json').write_text(json.dumps(parity_inputs,indent=2)+'\n')
    (pd/'SPEC.md').write_text((PARENT/'independent_parity_review_v1/SPEC.md').read_text()+'''

## Periodic source prototype binding

This is the exact existing parity comparator, now bound to the periodic-source
prototype's vacuum-only 57 results. The original two-attempt baseline and all
arithmetic/printing/derivative rules are unchanged. This adapter was prepared by
the prototype author. Independent actual numerical admission remains separate;
no review was executed during preparation and no periodic claim is supported.
''')
    assert (pd/'COMPARE.py').read_bytes()==compare.read_bytes()
    ready_common={'new_native_calls':0,'review_executed':False,'result':pin(result_path),
        'source_manifest':cfg['source_manifest'],'source_preparation':pin(__file__),
        'environment':{k:'1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS',
                    'MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','PYTHONDONTWRITEBYTECODE')}}
    for folder,source_name,argv,extra in (
        (ad,'REVIEW_ACTUAL.py',['--config',str(ad/'CONFIG_RUN.json'),'--out',
              str(WORK/'independent_actual_review_v1/actual_v1')],
              {'configuration':pin(ad/'CONFIG_RUN.json'),'numerical_algorithm_unchanged':True}),
        (pd,'COMPARE.py',['--result',str(result_path),'--out',
              str(WORK/'independent_parity_review_v1/actual_v1')],
              {'inputs':pin(pd/'INPUTS.json'),'source_byte_identical_to_parent':True})):
        ready=dict(ready_common,status='saved_review_prepared_unexecuted',
                   implementation=pin(folder/source_name),specification=pin(folder/'SPEC.md'),
                   argv=['/home/jwestrob/.pyenv/versions/3.12.0/bin/python3',str(folder/source_name)]+argv,
                   **extra)
        (folder/'READY.json').write_text(json.dumps(ready,indent=2)+'\n')
    report={'status':'prepared_existing_saved_reviewers_no_execution',
        'actual_ready':pin(ad/'READY.json'),'parity_ready':pin(pd/'READY.json'),
        'result':pin(result_path),'actual_algorithm_reversibly_identical':True,
        'parity_source_byte_identical':True,'new_native_calls':0,'new_numerical_replay':False,
        'cost':{'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    (HERE/'REVIEWERS_READY.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'ready':pin(HERE/'REVIEWERS_READY.json'),**report}))


if __name__=='__main__':main()
