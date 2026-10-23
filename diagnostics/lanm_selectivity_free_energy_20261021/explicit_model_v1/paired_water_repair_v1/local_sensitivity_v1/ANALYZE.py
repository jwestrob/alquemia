#!/usr/bin/env python3
"""Saved four-parameter Jacobian diagnosis. No molecular program or new candidate."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import time
import numpy as np

THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
KINDS, GEOMETRIES = ('La', 'Dy', 'paired'), ('water_central', 'aquo8_central')


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    if pin(record['path']) != record:
        raise ValueError('Changed evidence: '+record['path'])
    return Path(record['path'])


def close(a, b):
    if not np.allclose(a, b, atol=1e-10, rtol=1e-10):
        raise ValueError('Saved weighted/raw bookkeeping does not close')


def raw_vector(score):
    values = [score['balanced_exchange_error_kcal_mol']]
    for kind in KINDS:
        for q in GEOMETRIES:
            values.extend(score['per_fixture'][kind][q]['raw_body_components_kcal_mol_A'])
    return np.asarray(values, dtype=float)


def main(config_path, output):
    start, cpu = time.monotonic(), time.process_time()
    cfg = json.loads(Path(config_path).read_text())
    assert cfg['implementation'] == pin(__file__) and all(os.environ.get(k) == '1' for k in THREADS)
    verify(cfg['specification'])
    result = json.loads(verify(cfg['parent_result']).read_text())
    assert cfg['parent_result']['sha256'] == '524d093fb920304e254369097a96843408ab21597aa7898f5d1a3947c64fe6b6'
    assert result['status'] == 'finite_pilot_complete' and result['candidate_frozen']
    parent_config = json.loads(verify(result['config']).read_text())
    for key in ('implementation', 'specification', 'objective', 'search', 'evaluator', 'evaluator_config', 'root_release'):
        verify(parent_config[key])
    search = result['search']; records = search['evaluations']; best = records[search['best_key']]
    assert best['parameters'] == result['candidate_parameters'] and best['score'] == result['selected_score']
    z = np.asarray(best['normalized_parameters'], dtype=float)
    physical = np.asarray(best['parameters'], dtype=float)
    widths = np.asarray([p['upper']-p['lower'] for p in result['family']['parameters']], dtype=float)
    r, raw = np.asarray(best['score']['weighted_residuals'], dtype=float), raw_vector(best['score'])
    # Reconstruct exact original block/category/body/translation/rotation weights.
    scale, slices, local_slices = [.5], {'exchange': (0, 1)}, {}
    offset = 1
    for kind in KINDS:
        begin = offset
        for q in GEOMETRIES:
            response = best['score']['per_fixture'][kind][q]['model_body_response']
            cats = sorted(set(b['category'] for b in response))
            qbegin = offset
            for category in cats:
                bodies = [b for b in response if b['category'] == category]
                for body in bodies:
                    terms = 1 if category == 'ion' else 2
                    coefficient = .5*np.sqrt(1/(len(cats)*len(bodies)*terms*3))/2/np.sqrt(2)
                    scale.extend([coefficient]*(terms*3)); offset += terms*3
            local_slices[kind+'__'+q] = (qbegin, offset)
        slices[kind] = (begin, offset)
    scale = np.asarray(scale)
    close(raw*scale, r)
    columns, raw_columns, sensitivity_keys = [], [], []
    for j, (lo, hi) in enumerate(search['sensitivity']['parameter_intervals_used']):
        endpoints = []
        for value in (lo, hi):
            target = z.copy(); target[j] = value
            found = [(key, v) for key, v in records.items() if v['normalized_parameters'] == target.tolist()]
            assert len(found) == 1 and found[0][1]['status'] == 'complete'
            endpoints.append(found[0]); sensitivity_keys.append(found[0][0])
        a, b = endpoints[0][1]['score'], endpoints[1][1]['score']
        columns.append((np.asarray(b['weighted_residuals'])-np.asarray(a['weighted_residuals']))/(hi-lo))
        raw_columns.append((raw_vector(b)-raw_vector(a))/(hi-lo))
    J, Jraw = np.column_stack(columns), np.column_stack(raw_columns)
    close(Jraw*scale[:, None], J)
    delta, _, rank, singular = np.linalg.lstsq(J, -r, rcond=1e-6)
    close(singular, search['sensitivity']['singular_values'])
    lower, upper = np.maximum(-.025, -z), np.minimum(.025, 1-z)
    factor = 1.
    for value, lo, hi in zip(delta, lower, upper):
        if value > 0: factor = min(factor, hi/value)
        if value < 0: factor = min(factor, lo/value)
    scaled = delta*factor
    bounded, bounded_loss, feasible_patterns = None, None, 0
    for pattern in itertools.product((-1, 0, 1), repeat=4):
        fixed = [j for j, state in enumerate(pattern) if state]
        free = [j for j, state in enumerate(pattern) if not state]
        trial = np.zeros(4)
        for j in fixed: trial[j] = lower[j] if pattern[j] < 0 else upper[j]
        target = -(r+J[:, fixed]@trial[fixed])
        if free: trial[free] = np.linalg.lstsq(J[:, free], target, rcond=1e-6)[0]
        if np.any(trial < lower-1e-12) or np.any(trial > upper+1e-12): continue
        trial = np.clip(trial, lower, upper); feasible_patterns += 1
        loss = float(np.linalg.norm(r+J@trial)**2)
        if bounded_loss is None or loss < bounded_loss:
            bounded, bounded_loss = trial, loss
    assert bounded is not None

    def summarize(weighted, raw_values):
        close(raw_values*scale, weighted)
        controls, per_fixture = {}, {}
        for kind in KINDS:
            a, b = slices[kind]
            rms, maximum = float(4*np.linalg.norm(weighted[a:b])), float(np.max(np.abs(raw_values[a:b])))
            controls[kind] = {'rigid_RMS_kcal_mol_A': rms, 'maximum_component_kcal_mol_A': maximum,
                              'predicted_pass': rms <= 2 and maximum <= 5}
        for label, (a, b) in local_slices.items():
            per_fixture[label] = {'rigid_RMS_kcal_mol_A': float(4*np.sqrt(2)*np.linalg.norm(weighted[a:b])),
                                  'maximum_component_kcal_mol_A': float(np.max(np.abs(raw_values[a:b])))}
        energy = float(raw_values[0])
        return {'objective': float(weighted@weighted),
            'weighted_block_losses': {k: float(weighted[a:b]@weighted[a:b]) for k, (a, b) in slices.items()},
            'balanced_energy_error_kcal_mol': energy, 'controls': controls, 'per_fixture': per_fixture,
            'all_predicted_controls_pass': abs(energy) <= 1 and all(v['predicted_pass'] for v in controls.values())}

    def prediction(d):
        rr, rrraw = r+J@d, raw+Jraw@d
        return {'normalized_delta': d.tolist(), 'physical_parameter_delta': (d*widths).tolist(),
            'hypothetical_parameters_not_evaluated': (physical+d*widths).tolist(),
            'max_normalized_component': float(np.max(np.abs(d))),
            'inside_original_parameter_box': bool(np.all(z+d >= 0) and np.all(z+d <= 1)),
            'summary': summarize(rr, rrraw), 'weighted_residual': rr.tolist()}
    comparisons = []
    for key, record in records.items():
        actual = record['score']; d = np.asarray(record['normalized_parameters'])-z
        pred = summarize(r+J@d, raw+Jraw@d)
        distance = float(np.max(np.abs(d)))
        comparisons.append({'key': key, 'role': record['role'], 'normalized_delta': d.tolist(),
            'max_normalized_component': distance, 'within_declared_local_neighborhood': distance <= .025,
            'used_for_J': key in sensitivity_keys, 'predicted': pred,
            'actual': {'objective': actual['objective'], 'balanced_energy_error_kcal_mol': actual['balanced_exchange_error_kcal_mol'],
                       'controls': actual['controls'], 'development_checks_passed': actual['development_checks_passed']},
            'actual_minus_predicted_objective': actual['objective']-pred['objective'],
            'actual_objective_below_frozen_candidate': actual['objective'] < best['score']['objective']})
    remainder = r+J@delta
    answer = {'status': 'saved_local_diagnosis_complete', 'implementation': pin(__file__), 'config': pin(config_path),
        'parent_result': cfg['parent_result'], 'selected_candidate_unchanged': best['parameters'],
        'baseline': summarize(r, raw), 'Jacobian_normalized': J.tolist(), 'sensitivity_keys': sensitivity_keys,
        'singular_values': singular.tolist(), 'rank': int(rank),
        'local_projection': {'orthogonal_remainder_weighted_loss': float(remainder@remainder),
            'within_span_weighted_loss': float((r-remainder)@(r-remainder)),
            'fraction_loss_outside_local_span': float((remainder@remainder)/(r@r)),
            'J_transpose_remainder_norm': float(np.linalg.norm(J.T@remainder)),
            'not_a_global_family_lower_bound': True},
        'unconstrained_Gauss_Newton': prediction(delta), 'uniformly_scaled_Gauss_Newton': prediction(scaled),
        'bounded_linear_least_squares': prediction(bounded), 'feasible_active_set_patterns': feasible_patterns,
        'comparison_to_all_saved_vectors': comparisons, 'new_molecular_calls': 0,
        'candidate_reselected': False, 'validation_outcomes_opened': False, 'physical_accuracy_qualified': False,
        'numpy_version': np.__version__, 'local_numerical_threads': 1,
        'elapsed_seconds': time.monotonic()-start, 'parent_CPU_seconds': time.process_time()-cpu}
    out = Path(output); out.mkdir(parents=True, exist_ok=False)
    (out/'RESULT.json').write_text(json.dumps(answer, indent=2, allow_nan=False)+'\n')
    lines = ['# Saved local parameter-response diagnosis', '',
        f"Selected candidate stays unchanged: `{best['parameters']}`. Zero new molecular calls.", '',
        f"Jacobian rank {rank}/4. Local orthogonal residual retains {answer['local_projection']['fraction_loss_outside_local_span']:.6f} of the current weighted loss; this is not a nonlinear-family bound.", '',
        '| Prediction | Loss | La maximum | Dy maximum | Paired maximum | Energy error | All predicted controls |',
        '|---|---:|---:|---:|---:|---:|---|']
    for name in ('unconstrained_Gauss_Newton', 'uniformly_scaled_Gauss_Newton', 'bounded_linear_least_squares'):
        pred = answer[name]['summary']; c = pred['controls']
        lines.append(f"| {name} | {pred['objective']:.6f} | {c['La']['maximum_component_kcal_mol_A']:.6f} | {c['Dy']['maximum_component_kcal_mol_A']:.6f} | {c['paired']['maximum_component_kcal_mol_A']:.6f} | {pred['balanced_energy_error_kcal_mol']:.6f} | {pred['all_predicted_controls_pass']} |")
    lines += ['', 'Maximum components have kcal/mol/Å units; energy errors are kcal/mol. The two bounded predictions satisfy max|Δz|≤0.025 and original bounds.', '',
        '| Existing local point | Role | Actual loss | Linear loss | Actual − linear |', '|---|---|---:|---:|---:|']
    for row in comparisons:
        if row['within_declared_local_neighborhood'] and row['key'] != search['best_key']:
            lines.append(f"| {row['key'][:12]} | {row['role']} | {row['actual']['objective']:.6f} | {row['predicted']['objective']:.6f} | {row['actual_minus_predicted_objective']:.6f} |")
    lines += ['', 'These saved one-coordinate samples test local linear accuracy but do not establish the actual behavior of a coupled step. Full directions, individual RMS controls, per-fixture predictions and all108 saved-point comparisons are in RESULT.json. No new point is evaluated or selected; root must declare any finite native test separately.', '']
    (out/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({'result': pin(out/'RESULT.json'), 'report': pin(out/'REPORT.md')}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--config', required=True); p.add_argument('--output', required=True)
    a = p.parse_args(); main(a.config, a.output)
