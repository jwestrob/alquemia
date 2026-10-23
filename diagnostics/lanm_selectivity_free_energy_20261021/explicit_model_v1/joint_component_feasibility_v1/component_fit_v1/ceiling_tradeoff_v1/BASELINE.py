#!/usr/bin/env python3
"""Saved auxiliary-epigraph certificate; no solver construction or solve."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time

import numpy as np


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def read(p):
    require(pin(p['path']) == p, 'Changed authority')
    return json.loads(Path(p['path']).read_text())


def module(p, name):
    require(pin(p['path']) == p, 'Changed source')
    spec = importlib.util.spec_from_file_location(name, p['path'])
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def load(inputs):
    base = module(inputs['comparison_implementation'], 'original_comparison')
    old = module(inputs['row_removal_implementation'], 'original_row_removal')
    require(pin(inputs['declaration']['path']) == inputs['declaration'], 'Exact declaration')
    full = read(inputs['full_result']); review = read(inputs['full_actual_review'])
    require(review['status'] == 'passed_saved_conic_actual_review'
            and review['certificate_checks_passed'] is True
            and review['result'] == inputs['full_result']
            and review['implementation'] == inputs['full_actual_reviewer'], 'Exact arithmetic admission')
    require(pin(inputs['full_actual_reviewer']['path']) == inputs['full_actual_reviewer'], 'Reviewer source')
    require(full['status'] == 'completed_saved_conic_comparison' and full['error'] is None
            and full['solver_calls'] == 1 and full['implementation'] == inputs['comparison_implementation'], 'Completed comparison only')
    require(review['solver_status'] == full['solver_status']
            and review['solver_success'] == full['solver_success'], 'Preserve actual solver status')
    cfg = read(full['configuration'])
    require(review['source_review'] == cfg['source_review'], 'Exact original source review')
    original_inputs, _, settings = base.context(cfg)
    require(full['origin_key'] == original_inputs['origin_key']
            and full['derivative_keys'] == original_inputs['derivative_keys'], 'Exact origin and24 probes')
    require(full['effective_settings'] == base.attributes(settings), 'Same effective settings')
    core = base.module(original_inputs['core'])
    p = old.problem_arrays(read(full['linearized_problem']))
    model = base.assemble(p)
    require(pin(full['conic_data']['path']) == full['conic_data'], 'Exact sparse matrix pin')
    with np.load(full['conic_data']['path'], allow_pickle=False) as a:
        for key, value in [('A_data', model['A'].data), ('A_indices', model['A'].indices),
                           ('A_indptr', model['A'].indptr), ('b', model['b']), ('c', model['c'])]:
            require(np.array_equal(a[key], value), 'Exact original sparse array '+key)
        require(tuple(a['A_shape']) == model['A'].shape, 'Exact matrix shape')
    require(read(full['row_layout']) == model['layout'], 'Exact row layout')
    return base, old, core, p, model, full, settings


def auxiliary_q(ratios):
    require(len(ratios) == 2 and all(math.isfinite(x) and x >= 0 for x in ratios), 'Two finite squared ratios')
    rho = max(ratios)
    cushion = 1e-8 * max(1., abs(rho))
    q = rho + cushion
    require(math.isfinite(q) and q > rho, 'Positive finite auxiliary cushion')
    return rho, cushion, q


def certify(p, full, core):
    x = np.asarray(full['solution']['x'], dtype=float)
    require(x.shape == (13,) and np.all(np.isfinite(x)) and x[-1] >= 0, 'Finite original vector')
    raw = core.certify_fraction(p, np.r_[x[:12], x[-1]**2], 1., full['solver_success'])
    ratios = [raw['certificate']['blocks'][k]['ratio'] for k in ('P_energy', 'B_energy')]
    rho, cushion, q = auxiliary_q(ratios)
    raised = core.certify_fraction(p, np.r_[x[:12], q], 1., False)
    c = raised['certificate']
    require(raw['parameters'] == raised['parameters'] and raw['key'] == raised['key'], 'Physical vector unchanged')
    require(c['in_box_and_trust'] and c['all25_predicted'] and c['fixed_old_ceiling_pass']
            and c['all_residual_margins_strictly_positive'] and len(c['residual_constraints']) == 851
            and raised['eligible'], 'Strict unchanged full-constraint baseline')
    return {'physical_direction': x[:12].tolist(), 'raw_solver_status': full['solver_status'],
            'raw_q': full['returned_q'], 'raw_certificate': raw,
            'actual_component_rho': rho, 'auxiliary_cushion': cushion, 'auxiliary_q': q,
            'certificate': raised, 'physical_vector_changed': False,
            'interpretation': 'Saved linearized feasibility only; no solver convergence or native admission.'}


def run(inputs_path, out):
    wall, cpu = time.monotonic(), time.process_time()
    inputs = json.loads(Path(inputs_path).read_text())
    require(inputs['baseline_implementation'] == pin(__file__), 'Exact baseline source')
    out = Path(out).resolve()
    require(str(out) == inputs['baseline_output_directory'], 'Fixed fresh baseline output')
    out.mkdir(parents=True, exist_ok=False)
    result = {'status': 'started', 'implementation': pin(__file__), 'inputs': pin(inputs_path),
              'solver_calls': 0, 'native_calls': 0, 'error': None}
    try:
        base, old, core, p, model, full, settings = load(inputs)
        result['full_result'] = inputs['full_result']; result['full_actual_review'] = inputs['full_actual_review']
        result['linearized_problem'] = full['linearized_problem']; result['conic_data'] = full['conic_data']
        result['baseline'] = certify(p, full, core)
        result['status'] = 'passed_saved_auxiliary_epigraph_baseline'
    except Exception as exc:
        result.update(status='failed_saved_auxiliary_epigraph_baseline', error=repr(exc))
    result['cost'] = {'wall_seconds': time.monotonic()-wall, 'CPU_seconds': time.process_time()-cpu, 'threads': 1}
    (out/'RESULT.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json'), 'error': result['error']}))
    return result


if __name__ == '__main__':
    a = argparse.ArgumentParser(); a.add_argument('--inputs', required=True); a.add_argument('--out', required=True)
    args = a.parse_args(); raise SystemExit(bool(run(args.inputs, args.out)['error']))
