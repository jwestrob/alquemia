#!/usr/bin/env python3
"""Saved-only review of the exact four epsilon probes; never invoke an engine."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_epsilon_diagnostic_v1'
RUN = BASE / 'execution_v1/run_v1'
OUT = BASE / 'execution_actual_review_v1/replay_v2'
checks = []
def check(ok, name):
    checks.append({'name': name, 'passed': bool(ok)})
    if not ok:
        raise AssertionError(name)
def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
def verify(value):
    check(pin(value['path']) == value, 'pin:' + value['path'])
    return Path(value['path'])
def load(value):
    return json.loads(verify(value).read_text())
def module(value, name):
    spec = importlib.util.spec_from_file_location(name, verify(value))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result
def no_native(*args, **kwargs):
    raise AssertionError('Native execution forbidden in this saved review')

start, cpu = time.monotonic(), time.process_time()
check(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')), 'one-thread controls')
OUT.mkdir(parents=True, exist_ok=False)
rp = pin(RUN / 'execution/RESULT.json')
r = load(rp)
check(r['status'] == 'epsilon_diagnostic_complete' and r['error'] is None and r['job_id'] == '17596', 'exact terminal job admitted')
cfg = r['config']
check(cfg == json.loads((RUN / 'CONFIG_RUN.json').read_text()), 'executed config exact')
check(cfg['implementation']['sha256'] == '416b7f5164ca813a08a5da5b4c5b512c77689c9c7ca39cc8a57e634014b6724b', 'frozen implementation')
impl = module(cfg['implementation'], 'epsilon_review_source')
manifest = load(r['manifest'])
prep = load(r['preparation_result'])
check(prep['manifest'] == r['manifest'] and prep['new_native_calls'] == 0, 'exact zero-call preparation')
release = load(cfg['root_release'])
check(release['status'] == 'approved_frozen_epsilon_diagnostic' and release['maximum_new_native_programs'] == 32 and release['workers'] == 4, 'finite root release')
adapter = module(cfg['authority_helper'], 'epsilon_review_authority')
data = load(cfg['original32'])
training = load(cfg['training'])
body = module(cfg['body_helper'], 'epsilon_review_body')
objective = module(cfg['objective'], 'epsilon_review_objective')
native = module(data['config']['native_helper'], 'epsilon_review_native')
native.invoke = no_native
receipt_count = fresh_count = reused_count = 0
maximum_component_closure = 0.
for probe in manifest['probes']:
    record = r['records'][probe['id']]
    check(record['status'] == 'complete' and set(record['model']) == set(impl.CELLS) and set(record['cells']) == set(impl.CELLS), probe['id'] + ':complete denominator')
    check(record['eta'] == probe['eta'] and record['epsilon'] == probe['epsilon'], probe['id'] + ':exact displacement')
    for gid, value in probe['reused_cells'].items():
        check(record['model'][gid] == value and record['cells'][gid]['reused'] and record['cells'][gid]['authority'] == value['provenance']['native_authority'], probe['id'] + ':' + gid + ':opposite-metal exact reuse')
        verify(value['provenance']['native_authority'])
        reused_count += 1
    for gid, inp in probe['inputs'].items():
        directory = RUN / 'execution' / probe['id'] / gid
        cell = record['cells'][gid]
        check(cell == json.loads((directory / 'RESULT.json').read_text()) and cell['status'] == 'complete' and cell['input'] == inp, probe['id'] + ':' + gid + ':cell receipt exact')
        for key, filename in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
            check((directory / filename).read_bytes() == verify(inp['files'][key]).read_bytes(), probe['id'] + ':' + gid + ':runtime input ' + key)
        check(len(cell['receipts']) == 2, probe['id'] + ':' + gid + ':two programs')
        for receipt, program, args in zip(cell['receipts'], ('analyze','testgrad'), (['system.xyz','E'],['system.xyz','Y','N'])):
            check(receipt['returncode'] == 0 and receipt['argv'] == [data['config']['executables'][program]['path']] + args and set(receipt['threads'].values()) == {'1'}, probe['id'] + ':' + gid + ':' + program + ':execution')
            adapter.program_admission(verify(receipt['output']), program, native)
            receipt_count += 1
        parsed = native.parse(directory / 'analyze.out', directory / 'testgrad.out', dict(training['rows'][gid]['input'], derivative_test=False))
        expected = {'energy': parsed['total_energy_kcal_mol'], 'forces': parsed['force_kcal_mol_A'], 'native': parsed}
        check(expected == cell['model'] == record['model'][gid], probe['id'] + ':' + gid + ':raw energy-force replay exact')
        check(parsed['finite_difference'] is None and parsed['polarization']['no_native_nonconvergence_warning'], probe['id'] + ':' + gid + ':native controls')
        maximum_component_closure = max(maximum_component_closure, abs(parsed['printed_component_sum_remainder_kcal_mol']))
        fresh_count += 1
    check(objective.evaluate(training, record['model'], body) == record['score'], probe['id'] + ':all objective entries replay exactly')
check((fresh_count,reused_count,receipt_count,r['native_program_invocations']) == (16,16,32,32), 'exact fresh/reused/program totals')
replayed = impl.analyze(manifest, r['records'], data, training, body, objective)
numeric_differences = []
def compare_saved(a, b, path='analysis'):
    if isinstance(a, dict):
        check(isinstance(b,dict) and set(a)==set(b), path + ':keys')
        for k in a: compare_saved(a[k],b[k],path+':'+k)
    elif isinstance(a,list):
        check(isinstance(b,list) and len(a)==len(b), path + ':length')
        for i,(x,y) in enumerate(zip(a,b)): compare_saved(x,y,path+':'+str(i))
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        delta=abs(a-b);numeric_differences.append(delta)
        check(delta <= 1e-12 + 1e-12*abs(b),path+':saved arithmetic replay precision')
    else: check(a==b,path+':exact')
compare_saved(replayed,r['analysis'])
a = r['analysis']
base = np.array(manifest['baseline']['score']['weighted_residuals'])
J4 = np.array(manifest['four_column_Jacobian'])
Je = np.column_stack([(np.array(r['records'][m+'_plus']['score']['weighted_residuals']) - np.array(r['records'][m+'_minus']['score']['weighted_residuals'])) / (r['records'][m+'_plus']['eta']-r['records'][m+'_minus']['eta']) for m in ('La','Dy')])
U, s, _ = np.linalg.svd(J4, full_matrices=False)
threshold = s[0] * 1e-6
B = U[:,s > threshold]
W = Je - B @ (B.T @ Je)
C, snew, _ = np.linalg.svd(W, full_matrices=False)
C = C[:,snew > threshold]
old_res = base - B @ (B.T @ base)
new_res = old_res - C @ (C.T @ base)
check(np.array_equal(Je, np.array(a['epsilon_columns_per_eta'])), 'independent epsilon columns')
check(abs(threshold-a['absolute_singular_threshold']) < 1e-14 and B.shape[1] == 4 and C.shape[1] == 2, 'same absolute threshold and actual increment rank')
check(np.max(np.abs(new_res - np.array(a['residual_components']['after_nested_six_span']))) < 1e-14, 'independent nested projection')
check(all(v['passed'] for v in a['linearity'].values()), 'all eight midpoint gates pass')
non_vdw = [abs(v) for item in a['linearity'].values() for name, entry in item['native_energy_components'].items() if name != 'Van der Waals' for v in entry.values()]
check(max(non_vdw) == 0., 'all non-vdW energy components unchanged including polarization')
for p in manifest['protected_parents']:
    verify(p)
verify(rp)
check(not r['new_candidate_selected'] and not r['reserved_validation_opened'] and not a['parameter_solution_computed'] and not a['physical_accuracy_qualified'], 'diagnostic-only scientific scope')
event = json.loads((RUN / 'completion_event.json').read_text())['records']['17596']
check(event['state'] == 'COMPLETED' and event['exit_code'] == '0:0' and int(event['allocated_cpus']) == 4, 'actual terminal scheduler cost')
loss = a['weighted_residual_losses']
receipt = {
    'status':'passed_actual_native_and_projection_review', 'implementation':pin(__file__), 'result':rp,
    'reviewed_source':cfg['implementation'], 'executed_config':pin(RUN / 'CONFIG_RUN.json'), 'preparation':r['preparation_result'], 'manifest':r['manifest'],
    'counts':{'fresh_cells':fresh_count,'opposite_metal_reused_rows':reused_count,'native_programs':receipt_count,'all_residual_entries':len(base)},
    'check_count':len(checks),'checks':checks,'original_singular_values':s.tolist(),'absolute_singular_threshold':float(threshold),'increment_singular_values':snew.tolist(),
    'maximum_saved_analysis_replay_absolute_difference':max(numeric_differences),'saved_arithmetic_precision_only':{'absolute':1e-12,'relative':1e-12,'scientific_criteria_unchanged':True},
    'weighted_residual_losses':loss,'additional_projected_reduction_fraction':(loss['after_original_four_span']-loss['after_nested_six_span'])/loss['after_original_four_span'],
    'weighted_residual_decomposition':a['weighted_residual_decomposition'], 'epsilon_column_norms':a['column_norms'],
    'maximum_midpoint_energy_error_kcal_mol':max(abs(v['midpoint_energy_error_kcal_mol']) for v in a['linearity'].values()),
    'maximum_midpoint_force_error_kcal_mol_A':max(v['maximum_midpoint_force_component_error_kcal_mol_A'] for v in a['linearity'].values()),
    'maximum_component_closure_kcal_mol':maximum_component_closure,'maximum_non_vdw_component_change':max(non_vdw),
    'scientific_interpretation':['Both epsilon columns add resolved directions outside the four-parameter local span at the frozen point.', 'The unconstrained local projection reduces weighted loss by about two thirds; no parameter solution or nonlinear fitted model was evaluated.', 'Projected aquo8 losses increase while most dimer losses fall; no individual control, transfer or physical-accuracy qualification follows.', 'This supports a separately declared bounded six-parameter development test, with all existing per-fixture controls retained.'],
    'native_cost':{'scheduler_elapsed_seconds':int(event['elapsed_seconds']),'allocated_CPUs':int(event['allocated_cpus']),'allocated_CPU_seconds':int(event['elapsed_seconds'])*int(event['allocated_cpus']),'runner_wall_seconds':r['elapsed_seconds'],'parent_CPU_seconds':r['parent_cpu_seconds'],'aggregate_child_user_CPU_seconds':r['aggregate_child_user_cpu_seconds'],'aggregate_child_system_CPU_seconds':r['aggregate_child_system_cpu_seconds'],'per_receipt_overlapping_CPU_counters_not_summed':True},
    'review_cost':{'wall_seconds':time.monotonic()-start,'process_CPU_seconds':time.process_time()-cpu,'maximum_threads':1,'new_native_calls':0},
    'new_candidate_selected':False,'new_fit_performed':False,'reserved_validation_opened':False
}
(OUT / 'REVIEW.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
print(json.dumps({'review':pin(OUT/'REVIEW.json'),'status':receipt['status'],'checks':len(checks),'losses':loss,'cost':receipt['review_cost']}))
