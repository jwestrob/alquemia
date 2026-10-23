#!/usr/bin/env python3
"""Saved-only continuation delta review; no engine imports or native calls."""
import ast
import concurrent.futures
import copy
import hashlib
import json
import os
from pathlib import Path
import threading
import time
import types
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = HERE.parent / 'paired_water_radial_continuation_v1'
WORK = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'
OUT = WORK / HERE.name
checks = []

def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def check(ok, name):
    if not ok:
        raise AssertionError(name)
    checks.append({'case': name, 'passed': True})

def read(p):
    check(pin(p['path']) == p, 'pin ' + p['path'])
    return json.loads(Path(p['path']).read_text())

def token(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()

def require(ok, message):
    if not ok:
        raise ValueError(message)

def batch_controls(driver):
    """Execute the exact extracted nested batch against fake E/F callbacks."""
    tree = ast.parse(driver.read_text())
    batch = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == 'batch')
    code = compile(ast.fix_missing_locations(ast.Module(body=[batch], type_ignores=[])), str(driver), 'exec')
    cells = tuple('cell' + str(i) for i in range(8))
    snapshots = []
    records = {}

    class FakeEvaluator:
        def __init__(self, failure=False):
            self.stop_event = threading.Event()
            self.requests = []
            self.calls = []
            self.failure = failure
            self.sibling_started = threading.Event()
            self.failure_saved = threading.Event()
        def evaluate(self, values):
            self.calls.append(values)
            request = {'vector': values, 'cells': {g: {'status': 'not_attempted'} for g in cells}}
            self.requests.append(request)
            if self.failure and values[0] == 1:
                require(self.sibling_started.wait(5), 'Fake sibling did not start')
                request['cells'][cells[0]] = {'status': 'complete', 'fake_authority': 'retained'}
                raise ValueError('deliberate saved-only failure')
            if self.failure:
                self.sibling_started.set()
                require(self.failure_saved.wait(5), 'Immediate failure snapshot was not made')
            request['cells'] = {g: {'status': 'complete'} for g in cells}
            return {'values': values}
        def snapshot(self):
            return {'requests': copy.deepcopy(self.requests)}

    class SmallExecutor(concurrent.futures.ThreadPoolExecutor):
        def __init__(self, max_workers):
            super().__init__(max_workers=min(2, max_workers))

    fake_concurrent = types.SimpleNamespace(futures=types.SimpleNamespace(
        ThreadPoolExecutor=SmallExecutor, as_completed=concurrent.futures.as_completed))
    ev = FakeEvaluator()
    def persist():
        snapshots.append(copy.deepcopy(records))
        if any(v['status'] == 'unavailable' for v in records.values()):
            ev.failure_saved.set()
    namespace = {'token': token, 'records': records, 'require': require, 'np': np,
        'lower': np.zeros(4), 'width': np.ones(4), 'CELLS': cells,
        'persist': persist, 'ev': ev, 'obj': types.SimpleNamespace(evaluate=lambda data, model, body: {'objective': model['values'][0]}),
        'body': None, 'concurrent': fake_concurrent}
    ev.data = None
    exec(code, namespace)
    points = [{'parameters': [1., 0., 0., 0.]}, {'parameters': [2., 0., 0., 0.]}]
    returned = namespace['batch'](points + points[:1], 'trial', 1)
    check(len(ev.calls) == 2 and len(returned) == 3, 'Duplicate full vector reuses the completed batch entry')
    check(len(snapshots) == 2, 'Successful batch writes full RESULT only at start and end')
    check(all(v['status'] == 'scheduled' and set(v['cells']) == set(cells) for v in snapshots[0].values()), 'Batch start declares every point and full eight-cell denominator')
    check(all(v['status'] == 'complete' for v in snapshots[-1].values()), 'Successful batch ends only after every unique point completes')

    records.clear(); snapshots.clear(); ev = FakeEvaluator(True); ev.data = None; namespace['ev'] = ev
    try:
        namespace['batch'](points, 'derivative', 2)
    except RuntimeError as exc:
        check('no next batch' in str(exc), 'Failed batch cannot return a selectable result')
    else:
        raise AssertionError('Failure not propagated')
    check(ev.stop_event.is_set(), 'First callback failure sets global stop')
    check(any(list(s.values())[0]['status'] == 'unavailable' and list(s.values())[1]['status'] == 'scheduled' for s in snapshots), 'Immediate failure write occurs before the in-flight sibling finishes')
    final = snapshots[-1]
    a, b = (final[token(p['parameters'])] for p in points)
    check(a['status'] == 'unavailable' and b['status'] == 'complete', 'Completed sibling retained after another worker fails')
    check(set(a['cells']) == set(cells) and a['cells'][cells[0]]['fake_authority'] == 'retained', 'Final failure snapshot retains partial cell progress and full denominator')
    check(len(snapshots) == 4, 'Failure path writes start, immediate failure, drained batch and detailed failure')
    return {'success_full_writes': 2, 'failure_full_writes': 4, 'maximum_fake_workers': 2, 'native_calls': 0}

def main():
    start = time.monotonic()
    check(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')), 'One numerical thread')
    cfg = json.loads((SOURCE / 'CONFIG_PREPARE.json').read_text())
    ecfg = read(cfg['evaluator_config'])
    prep_path = WORK / 'paired_water_radial_continuation_v1/preparation_v2/RESULT.json'
    prep = read(pin(prep_path))
    prior = read(cfg['prior_result'])
    prior_cfg = read(prior['config'])
    prior_ecfg = read(prior_cfg['evaluator_config'])
    for key in ('implementation','specification','evaluator','objective','body_helper','water48','center_objective','joint_math','training_result'):
        check(pin(cfg[key]['path']) == cfg[key], 'Current bound source ' + key)
    for key in ('objective','body_helper','water48','center_objective','joint_math','training_result'):
        check(cfg[key] == prior_cfg[key], 'Unchanged science dependency ' + key)
    for key in ('resolver','native_helper','executables','source_inputs','baseline_inputs','fixed_pair_epsilon','initial_vector'):
        check(ecfg[key] == prior_ecfg[key], 'Unchanged physical context ' + key)
    check(prep['status'] == 'prepared' and prep['error'] is None and prep['new_native_program_invocations'] == 0, 'Actual v2 preparation passed without native calls')
    check(prep['config'] == pin(SOURCE / 'CONFIG_PREPARE.json') and prep['implementation'] == cfg['implementation'], 'Preparation binds final driver/config')
    check(prep['initial']['parameters'] == prior['selected']['parameters'] and prep['initial']['score'] == prior['selected']['score'], 'Exact17589 selected start and unchanged objective replay')
    check(prep['initial_center_diagnostics'] == prior['selected_center_diagnostics'], 'Original-center controls and method sensitivity unchanged')
    for gid, model in prep['initial']['model'].items():
        old = prior['selected']['model'][gid]
        check(all(model[k] == old[k] for k in ('energy','forces','native')), 'Exact starting native E/F ' + gid)
        check(model['provenance']['native_authority'] == old['provenance']['native_authority'] and model['provenance']['native_implementation'] == ecfg['prior_evaluator'] and model['provenance']['inherited_native'], 'Original native authority/source retained ' + gid)
    check(prep['evaluator']['inherited_native_cells_indexed'] == 672 and prep['evaluator']['inherited_native_cells_admitted'] == 8, 'Complete old native inventory indexed; eight actual start records replayed')
    check(ecfg['maximum_parameter_vectors'] == 385 and ecfg['maximum_native_program_invocations'] == 4096 and ecfg['maximum_native_concurrency'] == 8, 'Additional independent work ceilings')
    check(prior['completed_outer_steps'] == 12 and prior['new_parameter_vectors'] == 132 and prior['new_native_program_invocations'] == 1344, 'Parent actual cost remains separate')
    for p in prep['protected_parents']:
        check(pin(p['path']) == p, 'Immutable protected parent ' + p['path'])
    controls_path = WORK / 'paired_water_radial_continuation_v1/saved_controls_v2/RESULT.json'
    controls = read(pin(controls_path))
    check(controls['status'] == 'passed' and controls['new_native_calls'] == 0 and controls['evaluator'] == cfg['evaluator'], 'Author cache/budget controls bind final evaluator with zero calls')
    check(all(x['passed'] for x in controls['records']), 'All narrow author cache controls passed')
    batch = batch_controls(Path(cfg['implementation']['path']))
    result = {'status': 'passed', 'reviewer': pin(__file__), 'sources': {k: cfg[k] for k in ('implementation','specification','evaluator','evaluator_config')},
        'preparation': pin(prep_path), 'author_controls': pin(controls_path), 'prior_result': cfg['prior_result'],
        'checks': checks, 'passed_checks': len(checks), 'batch_failure_simulation': batch,
        'new_native_calls': 0, 'maximum_local_workers': 2, 'numerical_threads': 1,
        'elapsed_seconds': time.monotonic() - start, 'reserved_validation_opened': False,
        'physical_accuracy_qualified': False, 'scope': 'start/budget/read-only reuse/runtime checkpoint delta only'}
    target = OUT / 'REVIEW_v2.json'
    require(not target.exists(), 'Fresh independent result')
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'checks': result['passed_checks'], 'result': pin(target), 'elapsed_seconds': result['elapsed_seconds']}))

if __name__ == '__main__':
    main()
