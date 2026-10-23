#!/usr/bin/env python3
"""Pure callback controls; never imports or invokes the molecular evaluator."""
import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import threading
import time


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(output):
    start, cpu = time.monotonic(), time.process_time()
    here = Path(__file__).resolve().parent
    pins = {}
    for name in ('OBJECTIVE.py', 'SEARCH.py', 'PILOT.py'):
        data = (here/name).read_bytes()
        ast.parse(data)
        pins[name] = hashlib.sha256(data).hexdigest()
    parameters = [{'initial': .5, 'lower': 0., 'upper': 1.} for _ in range(4)]
    target = [.3, .4, .6, .7]
    def quadratic(x):
        # Reorder callback completion deliberately; selection must use poll order.
        time.sleep(.0001*(1+int(x[0]*100)%3))
        residual = [a-b for a, b in zip(x, target)]
        return {'model': {}, 'score': {'objective': sum(v*v for v in residual),
                'weighted_residuals': residual, 'development_checks_passed': False}}
    search = load(here/'SEARCH.py', 'source_control_quadratic')
    state = search.search(parameters, quadratic, threading.Event(), lambda state: None, workers=4)
    assert state['status'] == 'finite_pilot_complete' and state['numerical_optimization_complete']
    assert all(abs(a-b) < 1e-10 for a, b in zip(state['candidate_parameters'], target))
    assert state['sensitivity']['rank'] == 4 and state['sensitivity']['candidate_not_reselected_from_sensitivity']
    for poll in state['polls']:
        assert poll['complete']
        ordered = [poll['origin']]+poll['trial_keys']
        expected = min(ordered, key=lambda k: state['evaluations'][k]['score']['objective'])
        assert poll['selected'] == expected
    cap = load(here/'SEARCH.py', 'source_control_small_cap')
    cap.SEARCH_CAP, cap.VECTOR_CAP = 9, 17  # In-memory synthetic check only.
    limited = cap.search(parameters, quadratic, threading.Event(), lambda state: None, workers=4)
    assert limited['optimizer_budget_limited'] and len(limited['polls']) == 1
    assert limited['search_distinct_vectors'] == 9 and limited['distinct_vectors'] <= 17
    assert limited['unstarted_poll']['reason'] == 'complete_poll_exceeds_remaining_search_budget'
    event, observed = threading.Event(), []
    def fail(x):
        if x != [.5]*4:
            event.set()
            raise RuntimeError('Deliberate callback failure; no molecular program')
        return quadratic(x)
    try:
        search.search(parameters, fail, event, lambda state: observed.append(json.loads(json.dumps(state))), workers=4)
        raise AssertionError('Failure was not propagated')
    except RuntimeError as exc:
        assert 'Vector/native failure' in str(exc)
    assert event.is_set() and observed[-1]['status'] == 'failed_or_unavailable'
    assert not observed[-1]['polls'] and any(r['status'] == 'unavailable' for r in observed[-1]['evaluations'].values())
    result = {'status': 'passed', 'source_pins': pins, 'local_workers_maximum': 4, 'new_native_program_invocations': 0,
        'controls': {'quadratic_local_stop': True, 'completion_order_does_not_select': True,
                     'complete_poll_budget_stop': True, 'global_callback_failure_stop': True,
                     'sensitivity_not_reselected': True},
        'quadratic_vectors': state['distinct_vectors'], 'small_cap_vectors': limited['distinct_vectors'],
        'elapsed_seconds': time.monotonic()-start, 'parent_CPU_seconds': time.process_time()-cpu}
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with Path(output).open('x') as handle:
        handle.write(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', required=True)
    main(p.parse_args().output)
