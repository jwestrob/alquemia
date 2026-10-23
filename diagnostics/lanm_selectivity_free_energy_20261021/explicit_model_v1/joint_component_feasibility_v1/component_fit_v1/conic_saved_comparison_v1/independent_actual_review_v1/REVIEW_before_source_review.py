#!/usr/bin/env python3
"""Read-only independent arithmetic review of the single saved conic solve.

Does not import the producer, CORE, an optimizer, or a molecular evaluator.
Floating consistency is distinct from solver success and native qualification.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time

import numpy as np
from scipy import sparse

THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
           'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')
OLD = ('old_EF', 'pair_mu', 'C_mu', 'chi')
COMPONENTS = ('P_energy', 'B_energy')
CHECKS = 0


def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def read(ref):
    require(pin(ref['path']) == ref, 'Exact artifact pin: ' + ref['path'])
    return json.loads(Path(ref['path']).read_text())


def near(a, b, label):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    require(a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), label + ': finite shape')
    # Arithmetic replay tolerance only: never used to turn a negative margin positive.
    require(np.allclose(a, b, rtol=2e-12, atol=2e-12), label + ': arithmetic replay')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def scalar_report(p, d, q, success):
    """Original scalar norms and physical constraints, without conic matrices."""
    def affine(block):
        return [float(r) + math.fsum(float(a)*float(x) for a, x in zip(row, d))
                for r, row in zip(block['r'], block['J'])]
    margins, norms, gates = {}, {}, {}
    for name, block in p['blocks'].items():
        values = affine(block)
        sq = math.fsum(v*v for v in values)
        ratio = sq / block['denominator']
        margins['block/' + name] = (q if name in COMPONENTS else 1.) - ratio
        norms[name] = {'predicted_squared_norm': sq, 'fixed_denominator': block['denominator'], 'ratio': ratio}
    for name, value in zip(p['energy']['names'], affine(p['energy'])):
        margins['energy/upper/' + name] = 1-value
        margins['energy/lower/' + name] = 1+value
        gates['energy/' + name] = {'error_kcal_mol': value, 'passed': abs(value) <= 1}
    for name, block in p['forces'].items():
        values = affine(block)
        rms2 = math.fsum(float(w)*v*v for w, v in zip(block['weights'], values))
        maximum = max(abs(v) for v in values)
        margins['force/RMS/' + name] = 1-rms2/4
        for i, value in enumerate(values):
            margins['force/upper/' + name + '/' + str(i)] = 1-value/5
            margins['force/lower/' + name + '/' + str(i)] = 1+value/5
        gates['force/' + name] = {'rigid_RMS_kcal_mol_A': math.sqrt(rms2),
                                'maximum_body_component_error_kcal_mol_A': maximum,
                                'passed': rms2 <= 4 and maximum <= 5}
    params = [float(c+w*y) for c, w, y in zip(p['center'], p['width'], d)]
    box = [y-float(lo) for y, lo in zip(d, p['bounds_lower'][:12])]
    box += [float(hi)-y for y, hi in zip(d, p['bounds_upper'][:12])]
    valid_box = (min(box) >= 0 and 0 <= q <= p['q_current'] and
                 all(float(lo) <= v <= float(hi) for v, lo, hi in zip(params, p['lower'], p['upper'])))
    minimum = min(margins.values())
    all25 = all(gate['passed'] for gate in gates.values())
    ceilings = all(norms[k]['ratio'] <= 1 for k in OLD)
    eligible = valid_box and all25 and ceilings and q < p['q_current'] and (minimum >= 0 if success else minimum > 0)
    return {'margins': margins, 'blocks': norms, 'original25': gates, 'parameters': params,
            'eligible': eligible, 'valid_box': valid_box, 'minimum_box_margin': min(box),
            'minimum_margin': minimum, 'all25': all25, 'ceilings': ceilings}


def check_matrix(p, A, b, c, layout):
    """Reconstruct each row directly from its named physical quantity."""
    require(A.shape == (1817, 13) and b.shape == (1817,) and c.shape == (13,), 'Fixed dimensions')
    near(c, [0.]*12+[1.], 'Objective t')
    coverage = np.zeros(1817, dtype=int)
    expected_names = set(scalar_report(p, [0.]*12, p['q_current'], True)['margins'])
    expected_names |= {'bound/' + side + '/' + str(i) for side in ('lower', 'upper') for i in range(13)}
    require(len(layout) == len(expected_names) and {row['name'] for row in layout} == expected_names, 'Exact layout names')
    for row in layout:
        name, start, stop = row['name'], row['start'], row['stop']
        require(0 <= start < stop <= 1817, 'Row interval')
        coverage[start:stop] += 1
        if name.startswith('block/') or name.startswith('force/RMS/'):
            block = p['blocks'][name[6:]] if name.startswith('block/') else p['forces'][name[10:]]
            values, J = np.asarray(block['r']), np.asarray(block['J'])
            component = name in ('block/P_energy', 'block/B_energy')
            weights = (np.full(len(values), 1/math.sqrt(block['denominator'])) if name.startswith('block/')
                       else np.sqrt(block['weights'])/2)
            top = np.zeros(13); top[-1] = -1 if component else 0
            expected_A = np.vstack([top, np.column_stack([-weights[:, None]*J, np.zeros(len(values))])])
            expected_b = np.r_[0. if component else 1., weights*values]
            require(row['kind'] == 'SOC', 'SOC row kind')
        else:
            expected_A = np.zeros((1, 13))
            side = name.split('/')[1]
            if name.startswith('bound/'):
                i = int(name.split('/')[-1])
                lo = p['bounds_lower'][i] if i < 12 else 0.
                hi = p['bounds_upper'][i] if i < 12 else math.sqrt(p['q_current'])
                expected_A[0, i] = -1 if side == 'lower' else 1
                expected_b = [-lo if side == 'lower' else hi]
            else:
                if name.startswith('energy/'):
                    i = p['energy']['names'].index(name.split('/', 2)[2])
                    r, J, scale = p['energy']['r'][i], p['energy']['J'][i], 1.
                else:
                    key, index = name.split('/', 2)[2].rsplit('/', 1)
                    i = int(index); block = p['forces'][key]
                    r, J, scale = block['r'][i], block['J'][i], 5.
                sign = 1 if side == 'upper' else -1
                expected_A[0, :12] = sign*np.asarray(J)/scale
                expected_b = [1-sign*r/scale]
            require(row['kind'] == 'nonnegative' and stop-start == 1, 'Linear row kind')
        near(A[start:stop].toarray(), expected_A, name + ' A')
        near(b[start:stop], expected_b, name + ' b')
    require(np.all(coverage == 1), 'Layout covers every row once')


def margins(layout, vector):
    out = {}
    for row in layout:
        v = vector[row['start']:row['stop']]
        out[row['name']] = float(v[0]-np.linalg.norm(v[1:])) if row['kind'] == 'SOC' else float(v[0])
    return out


def run(result_path, output):
    tick, cpu = time.monotonic(), time.process_time()
    require(all(os.environ.get(k) == '1' for k in THREADS), 'One numerical thread')
    result_pin = pin(result_path); result = read(result_pin)
    require(result['status'] == 'completed_saved_conic_comparison' and result['error'] is None, 'Completed comparison')
    output = Path(output); output.mkdir(parents=True, exist_ok=False)
    report = {'status': 'started', 'implementation': pin(__file__), 'result': result_pin,
              'certificate_checks_passed': False, 'new_native_calls': 0, 'new_solver_calls': 0,
              'solver_success': False, 'error': None}
    try:
        cfg = read(result['configuration']); source = read(cfg['source_review'])
        report['source_review'] = cfg['source_review']
        require(read(source['configuration']) == read(pin(Path(__file__).resolve().parents[1]/'CONFIG_UNBOUND.json')), 'Reviewed configuration authority')
        require(pin(cfg['implementation']['path']) == cfg['implementation'] == result['implementation'] == source['implementation'], 'Producer source pins')
        fingerprint = digest({k: v for k, v in cfg.items() if k not in ('terminal_campaign', 'source_review', 'root_release')})
        require(source['status'] == 'passed_conic_saved_comparison_source_review' and
                fingerprint == source['configuration_fingerprint'] == result['configuration_fingerprint'], 'Source review binding')
        release = read(cfg['root_release']); fit = read(cfg['terminal_campaign']); inputs = read(cfg['inputs'])
        require(release['status'] == 'released_one_saved_origin24_conic_comparison' and
                release['terminal_campaign'] == cfg['terminal_campaign'] == result['terminal_campaign'] and
                release['source_review'] == cfg['source_review'] and release['configuration_fingerprint'] == fingerprint and
                release['maximum_solver_calls'] == 1 and release['new_native_calls'] == 0, 'Single comparison release')
        require(fit['status'] in ('bounded_component_fit_complete', 'failed_component_fit') and fit['mode'] == 'run' and 'cost' in fit, 'Terminal fit')
        require(result['origin_key'] == inputs['origin_key'] and result['derivative_keys'] == inputs['derivative_keys'], 'Frozen origin/probes')
        require(result['solver_calls'] == 1 and result['new_native_calls'] == 0 and result['new_candidate_evaluations'] == 0 and
                result['selected_candidate'] is None and result['automatic_fit_extension'] is False, 'Saved-only execution counts')
        require(result['effective_settings'] == read(cfg['settings']), 'Exact settings')
        marker = json.loads(Path(cfg['one_attempt_marker']).read_text())
        require(marker['configuration'] == result['configuration'] and marker['terminal_campaign'] == cfg['terminal_campaign'], 'Persistent single attempt marker')
        p, layout = read(result['linearized_problem']), read(result['row_layout'])
        require(pin(result['conic_data']['path']) == result['conic_data'], 'Matrix data pin')
        with np.load(result['conic_data']['path'], allow_pickle=False) as data:
            A = sparse.csc_matrix((data['A_data'], data['A_indices'], data['A_indptr']), shape=tuple(data['A_shape']))
            b, c = data['b'].copy(), data['c'].copy()
        require(tuple(p['blocks']) == OLD+COMPONENTS and len(p['forces']) == 15 and len(p['energy']['r']) == 10, 'Six blocks /25 gates')
        check_matrix(p, A, b, c, layout)
        x, s, z = (np.asarray(result['solution'][key], dtype=float) for key in ('x', 's', 'z'))
        require(x.shape == (13,) and s.shape == z.shape == (1817,) and all(np.all(np.isfinite(v)) for v in (x, s, z)), 'Finite solution')
        require(x[-1] >= 0, 'Nonnegative t before squaring')
        success = result['solver_status'] == 'Solved'
        require(result['solver_success'] == success and result['solution']['status'] == result['solver_status'], 'Truthful solver status')
        stored = result['independent_linear_residuals']
        primal, dual = A@x+s-b, A.T@z+c
        pcost, dcost = float(c@x), float(-b@z)
        values = {'primal_equation_residual': primal, 'dual_stationarity_residual': dual,
                  'primal_residual_inf': np.max(abs(primal)), 'dual_residual_inf': np.max(abs(dual)),
                  'primal_objective_t': pcost, 'dual_objective_t': dcost, 'signed_gap_t': pcost-dcost,
                  'complementarity_dot': float(s@z)}
        for key, value in values.items():
            near(value, stored[key], key)
        for key, vector in [('primal_affine_cone_margins', b-A@x), ('raw_slack_cone_margins', s), ('dual_cone_margins', z)]:
            replay = margins(layout, vector)
            require(set(replay) == set(stored[key]), key + ' names')
            for name, value in replay.items():
                near(value, stored[key][name], key + '/' + name)
        require(stored['rigorously_rounded_dual_certificate'] is False, 'Floating certificate label')
        near(result['returned_q'], x[-1]**2, 'Returned q')
        require([v['fraction'] for v in result['fractions']] == [1., .5, .25, .125], 'Four declared fractions')
        summaries = []
        for slot in result['fractions']:
            f = slot['fraction']; d = [f*float(v) for v in x[:12]]
            q = (1-f)*p['q_current']+f*float(x[-1]**2)
            replay = scalar_report(p, d, q, success); cert = slot['certificate']
            require(len(replay['margins']) == 851 and len(replay['original25']) == 25, 'Exact original constraints')
            for key in ('residual_constraints', 'blocks', 'original25'):
                expected = replay[{'residual_constraints': 'margins', 'blocks': 'blocks', 'original25': 'original25'}[key]]
                # Same explicitly ordered scalar arithmetic permits exact sign and gate replay.
                require(cert[key] == expected, 'Exact scalar certificate: ' + key)
            require(slot['parameters'] == replay['parameters'] and slot['eligible'] == replay['eligible'], 'Candidate identity / strict eligibility')
            require(cert['q_fraction'] == q and cert['in_box_and_trust'] == replay['valid_box'] and
                    cert['minimum_residual_constraint_margin'] == replay['minimum_margin'] and
                    cert['all25_predicted'] == replay['all25'] and cert['fixed_old_ceiling_pass'] == replay['ceilings'], 'Certificate summaries')
            require(slot['native_evaluated'] is False, 'No unevaluated candidate promotion')
            summaries.append({'fraction': f, 'eligible': replay['eligible'], 'q': q, 'minimum_margin': replay['minimum_margin'],
                              'all25_predicted': replay['all25'], 'fixed_ceilings': replay['ceilings']})
        report.update(status='passed_saved_conic_actual_review', certificate_checks_passed=True,
                      solver_status=result['solver_status'], solver_success=success,
                      returned_q=result['returned_q'], primal_residual_inf=float(np.max(abs(primal))),
                      dual_residual_inf=float(np.max(abs(dual))), signed_gap_t=pcost-dcost,
                      fractions=summaries, interpretation='Saved arithmetic and provenance verified. Solver success is separate; no rigorous optimality, nonlinear feasibility, native candidate admission or fit extension.')
    except Exception as exc:
        report.update(status='failed_saved_conic_actual_review', error=repr(exc))
    report['checks'] = CHECKS
    report['cost'] = {'wall_seconds': time.monotonic()-tick, 'CPU_seconds': time.process_time()-cpu, 'numerical_threads': 1}
    (output/'REVIEW.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': report['status'], 'review': pin(output/'REVIEW.json'), 'error': report['error']}))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--result', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    raise SystemExit(bool(run(args.result, args.out)['error']))
