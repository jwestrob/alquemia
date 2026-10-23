"""One saved-metadata review. No numpy/scipy, native calls, matrix loads or replay.

Scope declared before execution: pin the completed initializer/assembly/factor
chain, inspect saved state checks, factor flags and ten scalar probe records,
and reconcile parent-only scheduler cost. Reuse admitted source and initializer
checks. Never reconstruct h/H, PME, sparse matrices, factors or solve vectors.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import time

BASE = Path(__file__).resolve().parents[7]
PR = Path('explicit_model_v1/radial_drive_shape_v1/periodic_runtime_v1')
W = BASE / 'workspaces/lanm_selectivity_free_energy_20261021' / PR
RUN = W / 'central_certificate_standard_v1/run_v1'
OUT = W / 'central_certificate_standard_v1/independent_actual_review_v1/actual_v1'
EXPECTED = {
    'configuration': '42a25a50940edaba2dd5a7202d73f5cb4a089262a2a6304fa81ec3131b631dae',
    'release': '1623553d419966a8012569b0397f325cf9b44289bfa443634e6cb81da7eb533b',
    'assembly': '6282ffcbd18daede92ce78c1f8014992e0589e83cb921c3c3099a5ac2dad67cc',
    'factor': '507c53f53ef24efb84406d1f50c0a7d9299c44752bba2dbad117cc2858f7dd4e',
    'implementation': 'c1f46d4ed21bd83ef6688b62660e6d31cb0ce3740d2254b15934cabf400bb58c',
    'manifest': '427013d1a664ae8fc74a36fb339bf09af5e4df8c5ecb167ec7386a0f4f5f03da',
    'initializer': '99dcae9737ad4a9d158db6952ad8688870a40c807048f33e9059f554e0a80c03',
    'build': '5f9d6ebaf0f54991fd15b0cbd01ce4732eb785301ecc2dcaa16575406af282f4',
}

def pin(p):
    p = Path(p).resolve()
    h = hashlib.sha256()
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return {'path': str(p), 'sha256': h.hexdigest()}

def read(p):
    return json.loads(Path(p).read_text())

def main():
    wall, cpu = time.monotonic(), time.process_time()
    OUT.mkdir(parents=True, exist_ok=False)
    checks, verified = [], {}
    rec = {'status': 'started', 'review_implementation': pin(__file__),
           'scope': __doc__, 'native_calls': 0, 'solver_calls': 0,
           'factorizations': 0, 'matrix_array_reads': 0, 'numerical_replays': 0,
           'error': None}
    def check(name, ok):
        checks.append({'name': name, 'passed': bool(ok)})
        if not ok:
            raise ValueError(name)
    def verify(p):
        if p['path'] not in verified:
            check('pin ' + p['path'], pin(p['path']) == p)
            verified[p['path']] = p
        else:
            check('consistent pin ' + p['path'], verified[p['path']] == p)
        return p
    def load(p):
        return read(verify(p)['path'])
    try:
        check('one-thread environment', all(os.environ.get(k) == '1' for k in
              ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')))
        result_pin = pin(RUN / 'execution/RESULT.json')
        r = load(result_pin)
        for k, sha in EXPECTED.items():
            check('fixed ' + k, r[k]['sha256'] == sha)
            verify(r[k])
        cfg, rel, a, f, m, init = (load(r[k]) for k in
                                  ('configuration','release','assembly','factor','manifest','initializer'))
        sub = read(RUN / 'SUBMISSION.json')
        event = read(RUN / 'completion_event.json')
        sched = event['records']['20844']
        for k in ('configuration','release','manifest'):
            check('submission ' + k, sub[k] == r[k])
        verify(sub['wrapper'])
        check('terminal parent scheduler', sched['job_id'] == sub['job_id'] == '20844'
              and sched['state'] == 'COMPLETED' and sched['exit_code'] == '0:0'
              and int(sched['allocated_cpus']) == 12)
        (OUT / 'SCHEDULER_PARENT_RECORD.json').write_text(json.dumps(sched, indent=2)+'\n')
        check('one state release', rel['configuration'] == r['configuration']
              and rel['implementation'] == r['implementation']
              and rel['state_id'] == 'central' and rel['output'] == str((RUN/'execution').resolve())
              and rel['maximum_assemblies'] == rel['maximum_numeric_factors'] == 1
              and rel['maximum_solve_probes'] == 10)
        source_review = load(rel['source_review'])
        check('reviewed certificate helpers', source_review['certificate_sources'] == cfg['certificate_sources'])
        for p in list(cfg['certificate_sources'].values()) + list(cfg['authorities'].values()) + list(m['helpers'].values()):
            verify(p)
        for k in ('manifest','build','initializer'):
            check('configuration ' + k, cfg[k] == r[k])
        verify(cfg['nominal_contract']); verify(cfg['build_admission'])
        ia = load(rel['initializer_root_admission'])
        check('initializer admitted', ia['status'] == 'admitted_periodic_v2_central_initialization'
              and ia['result'] == r['initializer'])
        check('current initializer only', init['status'] == 'stage2_initialization_complete'
              and init['manifest'] == r['manifest'] and init['build'] == r['build']
              and init['counts'] == dict(native_programs=1,direct_initializations=1,
                  mutual_applications=0,induced_solves=0,gradient_calls=0,energy_only_calls=0))
        for k in ('expected_state','pme_export','candidate','source_manifest'):
            verify(r[k])
        check('state and PME exact outputs', r['expected_state'] == init['state_output']
              and r['pme_export'] == init['pme_export'])
        for k in ('xyz','key','parameters'):
            verify(r['inputs'][k])
        check('same literal identity through both phases', r['identity'] == init['identity'] == a['identity'] == f['identity'])
        check('cutoff12 V2 source identity', r['identity']['schema'] == 'PERIODIC_PARITY_STATE_V2'
              and r['identity']['settings']['ewaldcut_A'] == a['cutoff_A'] == 12.
              and r['identity']['cell'] == a['cell_A'] == [49.324]*3
              and r['candidate'] == m['model']['candidate'])
        check('no circular downstream admission', cfg['central_parity_admission'] is None
              and cfg['central_scout_review'] is None)
        check('successful supervisor', r['status'] == 'stage2_state_nominal_certificate_passed'
              and r['error'] is None and r['sufficient_numerical_test_passed']
              and r['native_calls'] == 0 and r['automatic_retry'] is False)
        for phase, obj in (('assembly',a),('factor',f)):
            for k in ('configuration','release','implementation','manifest','build','initializer','inputs','source_inputs','source_manifest'):
                check(phase+' join '+k, obj[k] == r[k])
            c = r['commands'][phase]
            check(phase+' one successful process', c['status'] == 'complete' and c['returncode'] == 0
                  and c['argv'][-2:] == ['--phase',phase] and c['PID'] > 0)
            for k in ('stdout','stderr'):
                verify(c[k]); check(phase+' empty '+k, Path(c[k]['path']).stat().st_size == 0)
            check(phase+' limits unchanged', obj['limits'] == cfg['caps'][phase])
            check(phase+' resource use', obj['cost']['wall_seconds'] < obj['limits']['wall']
                  and obj['cost']['CPU_seconds'] < obj['limits']['CPU']
                  and obj['cost']['peak_RSS_KiB']*1024 < obj['limits']['AS'])
            check(phase+' exact allocation', obj['allocation']['job_id'] == '20844'
                  and obj['allocation']['allocated_CPUs'] == 12 and obj['allocation']['qos'] == 'standard'
                  and obj['allocation']['memory_grant_bytes'] == 48000*1024**2
                  and not obj['allocation']['requeue'])
        check('one assembly, no hidden computation', a['status'] == 'stage2_assembly_complete'
              and a['error'] is None and a['assembly_attempts'] == 1
              and all(a[k] == 0 for k in ('native_calls','symbolic_attempts','numeric_attempts','solve_attempts')))
        check('full saved body checks', len(a['state_checks']) == 4 and all(
              c['passed'] and math.isfinite(c['max_abs_error']) and 0 <= c['max_abs_error'] <= c['bound']
              and c['bound'] == 1e-11*c['scale'] for c in a['state_checks']))
        check('exact sparse size and normalization', a['shape'] == [35874,35874]
              and a['active_site_count'] == 11958 and a['nnz'] == 38813535
              and a['info']['triangle'] == 'lower' and a['info']['normalization'] == 'sqrt(alpha) D sqrt(alpha)'
              and a['info']['mutual_scales_all_one'] and not a['info']['PME_included']
              and not a['info']['native_operator_parity_admitted'])
        check('factor input exact assembly', f['assembly'] == r['assembly'])
        norms = init['PME_checks']['operator_norms']
        check('admitted PME scopes', init['PME_checks']['site_checks']['passed']
              and init['PME_checks']['summary_comparison']['passed']
              and norms['q_negative_count'] == 0 and norms['q_min'] >= 0)
        n = f['operator']
        check('nominal shift exact scalar formula', n == f['nominal'] == r['nominal']
              and n['U_B'] == a['U_B'] == a['info']['unshifted_site_upper']
              and n['U_R'] == a['U_R'] == norms['U_R_nominal']
              and n['delta_odd'] == a['delta_odd'] == norms['evenness_perturbation_nominal']
              and n['U_nominal'] == n['U_B']+n['U_R']
              and n['tau'] == 1e-6*n['U_nominal']+n['delta_odd'])
        check('nominal-only scope', n['scope'] == 'nominal_saved_operator_only'
              and n['e_app'] is None and not n['application_error_bound_established'])
        check('one symbolic and numeric factor, ten solves', f['symbolic_attempts'] == f['numeric_attempts'] == 1
              and f['solve_attempts'] == 10 and f['assembly_attempts'] == f['native_calls'] == 0)
        check('factor pass without warnings', f['status'] == 'passed_nominal_static_sufficient_numerical_test'
              and f['error'] is None and f['backend_exception'] is None
              and not f['backend_warnings'] and not f['audit_warnings'])
        check('actual twelve BLAS and one OpenMP', f['actual_threads'] == {'openblas':12,'openmp':1})
        pre, post, audit = f['pre_numeric'],f['post_numeric'],f['audit']
        check('qualified backend state transition', pre['is_super'] and not pre['is_numeric']
              and post['is_numeric'] and post['is_ll'] and not post['is_super'] and post['is_lower']
              and audit['before_copy'] == audit['after_copy'] == post
              and all(x['raw_minor'] is None and x['raw_status'] is None for x in (pre,post)))
        names = [f'uniform_axis_{k}' for k in range(3)]+[f'metal_basis_{k}' for k in range(3)]+[f'rademacher_{k}' for k in range(4)]
        check('ten exact saved probe identities', [p['name'] for p in audit['probes']] == names)
        check('inherited arithmetic bound', audit['arithmetic_bound'] == 100*35874*2.220446049250313e-16)
        for p in audit['probes']:
            check('saved probe '+p['name'], p['finite'] and p['passed'] and all(
                  math.isfinite(p[k]) and 0 <= p[k] <= audit['arithmetic_bound']
                  for k in ('factor_action_error','solve_backward_error')))
        check('saved audit complete', audit['passed'] and audit['finite_positive_diagonal']
              and audit['original_factor_flags_preserved'] and audit['sparse_factor_shape'] == a['shape']
              and sorted(audit['permutation']) == list(range(35874))
              and audit['maximum_factor_action_error'] == max(p['factor_action_error'] for p in audit['probes'])
              and audit['maximum_solve_backward_error'] == max(p['solve_backward_error'] for p in audit['probes']))
        child_cpu = sum(c['cost']['child_user_CPU_seconds']+c['cost']['child_system_CPU_seconds'] for c in r['commands'].values())
        rec.update(status='passed_periodic_v2_central_certificate_actual_review',
                   result=result_pin, configuration=r['configuration'], release=r['release'],
                   implementation=r['implementation'], initializer=r['initializer'],
                   initializer_admission=rel['initializer_root_admission'], manifest=r['manifest'],
                   assembly=r['assembly'], factor=r['factor'], source_review=rel['source_review'],
                   scheduler_parent=pin(OUT/'SCHEDULER_PARENT_RECORD.json'),
                   counts={'assemblies':1,'symbolic_factors':1,'numeric_factors':1,'solve_probes':10,'native_calls':0},
                   nominal=n, state_checks=a['state_checks'], metal_stiffness=a['metal_stiffness'],
                   arrays_metadata_only=a['arrays'], array_payloads_not_recomputed=a['payload_sha256'],
                   backend={'pre':pre,'post':post,'warnings':[],'audit_warnings':[],
                            'factor_nnz':audit['sparse_factor_nnz'],'bound':audit['arithmetic_bound'],
                            'maximum_action_error':audit['maximum_factor_action_error'],
                            'maximum_solve_error':audit['maximum_solve_backward_error']},
                   actual_cost={'supervisor':r['cost'],'assembly':a['cost'],'factor':f['cost'],
                       'child_CPU_seconds_from_parent':child_cpu,
                       'total_measured_CPU_seconds':r['cost']['controller_CPU_seconds']+child_cpu,
                       'scheduler_parent_elapsed_seconds':int(sched['elapsed_seconds']),
                       'scheduler_allocated_CPU_seconds':int(sched['elapsed_seconds'])*int(sched['allocated_cpus'])},
                   conclusion='Passed the nominal sufficient numerical test for this exact saved cutoff12 radial14 central state only.',
                   limitations=['No sparse array/hash payload, h/H, PME, factor action or solve replay in this review; producer checks are inspected and pinned.',
                       'Sampled backward checks are not a whole-factor error norm, spectral forward-error bound or rigorous floating-application envelope.',
                       'Raw backend minor/status are unavailable. Symbolic supernodal to numeric simplicial LL storage is the qualified wrapper behavior.',
                       'Native full-operator parity, periodic energy/force consistency, physical accuracy and trajectory stability are not admitted by this result; e_app remains unknown and optional.'])
    except BaseException as exc:
        rec.update(status='failed_periodic_v2_central_certificate_actual_review', error=repr(exc))
    rec.update(checks=checks, checks_passed=sum(c['passed'] for c in checks), verified_pins=list(verified.values()),
               review_cost={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,
                            'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'threads':1})
    (OUT/'REVIEW.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps({'review':pin(OUT/'REVIEW.json'),'status':rec['status'],'error':rec['error'],
                      'checks_passed':rec['checks_passed'],'review_cost':rec['review_cost']}))
    return 0 if rec['error'] is None else 2

if __name__ == '__main__':
    raise SystemExit(main())
