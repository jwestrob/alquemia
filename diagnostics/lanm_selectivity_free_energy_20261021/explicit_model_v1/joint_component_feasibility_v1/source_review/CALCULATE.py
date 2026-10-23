#!/usr/bin/env python3
"""Read frozen vectors; diagnose saved component response. No evaluator/native API."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

import numpy as np
import scipy

HERE = Path(__file__).resolve().parent
THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
SLICES = {'old_EF': (0, 413), 'pair_mu': (413, 515), 'C_mu': (515, 527), 'chi': (527, 531)}
SCORE_KEYS = {'old_EF': 'old_EF', 'pair_mu': 'pair', 'C_mu': 'C', 'chi': 'chi'}
ERROR_ORDER = ('P8', 'B8', 'P9', 'B9')
PARENTS = tuple(f'aquo{n}_central_{m}' for n in (8, 9) for m in ('La', 'Dy'))


def require(ok, msg):
    if not ok:
        raise ValueError(msg)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def read(record):
    require(pin(record['path']) == record, 'Changed artifact: '+record['path'])
    return json.loads(Path(record['path']).read_text())


def token(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def close(a, b, label):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    require(a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), label+' shape/finite')
    require(np.all(np.abs(a-b) <= 1e-10*np.maximum(1., np.maximum(np.abs(a), np.abs(b)))), label+' arithmetic')
    return float(np.max(np.abs(a-b), initial=0.))


def save(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def gates(score):
    old = score['old_EF']
    result = {}
    for group, entries in old['energy_controls'].items():
        for name, row in entries.items():
            result['energy/'+group+'/'+name] = dict(row)
    for group, entries in old['per_fixture'].items():
        for name, row in entries.items():
            result['force/'+group+'/'+name] = {k: row[k] for k in ('rigid_RMS_kcal_mol_A', 'maximum_body_component_error_kcal_mol_A', 'passed')}
    require(len(result) == 25 and sum(k.startswith('energy/') for k in result) == 10, 'Original10+15 gate denominator')
    passed = all(row['passed'] is True for row in result.values())
    require(passed == score['original_EF_checks_passed'] == old['development_checks_passed'], 'Stored25 gate summary')
    return result


def cycles(decomposition, targets):
    group = targets['groups']['pair_and_collective_dipoles']
    require(set(group['parents']) == set(PARENTS) and len(group['pair_cycles']) == 34, 'Literal34/fourparent denominator')
    terms = {}
    for parent in PARENTS:
        n, metal = int(parent[4]), parent.rsplit('_', 1)[1]
        pairs = [r for r in group['pair_cycles'] if r['parent'] == parent]
        require({r['water_index'] for r in pairs} == set(range(1, n+1)) and len(pairs) == n, 'Literal parent water indices')
        C = {parent: 1, 'ion_'+metal: -1, parent+'__shell': -1}
        P = {}
        for row in pairs:
            for key, coefficient in row['model_terms'].items():
                P[key] = P.get(key, 0)+coefficient
        require(P['ion_'+metal] == -n, 'N-fold bare ion subtraction')
        B = dict(C)
        for key, coefficient in P.items():
            B[key] = B.get(key, 0)-coefficient
        terms[parent] = {'C': C, 'P': P, 'B': B}
        require(terms[parent] == decomposition['parents'][parent]['model_terms'], 'Admitted decomposition coefficients')
        require(decomposition['parents'][parent]['reference_electronic_plus_D4_kcal_mol'] is None, 'Absolute-parent target remains unavailable')
    refs = {}
    for n in (8, 9):
        row = decomposition['paired_Dy_minus_La']['W'+str(n)]
        refs[n] = row['reference_total_kcal_mol']
        for k in ('C', 'P', 'B'):
            close(refs[n][k], row['reference_electronic_kcal_mol'][k]+row['reference_D4_kcal_mol'][k], 'Separate electronic/D4 reference')
        close(refs[n]['C'], refs[n]['P']+refs[n]['B'], 'Reference C=P+B')
    return terms, refs


def component_values(energies, terms, refs):
    parents, paired = {}, {}
    closure = 0.
    for parent, channels in terms.items():
        row = {k: math.fsum(c*energies[x] for x, c in coeffs.items()) for k, coeffs in channels.items()}
        closure = max(closure, close(row['C'], row['P']+row['B'], 'Model C=P+B'))
        parents[parent] = row
    errors = []
    for n in (8, 9):
        la, dy = (parents[f'aquo{n}_central_{m}'] for m in ('La', 'Dy'))
        model = {k: dy[k]-la[k] for k in ('C', 'P', 'B')}
        err = {k: model[k]-refs[n][k] for k in model}
        close(err['C'], err['P']+err['B'], 'Paired error C=P+B')
        paired['W'+str(n)] = {'model_kcal_mol': model, 'error_kcal_mol': err}
        errors.extend((err['P'], err['B']))
    return parents, paired, np.asarray(errors), closure


def blocks(residual, component, jacobian=None, component_jacobian=None):
    result = {}
    for name, (a, b) in SLICES.items():
        result[name] = {'residual': residual[a:b]}
        if jacobian is not None:
            result[name]['jacobian'] = jacobian[a:b]
    for name, ix in (('P_energy', [0, 2]), ('B_energy', [1, 3])):
        result[name] = {'residual': component[ix]/math.sqrt(2.)}
        if component_jacobian is not None:
            result[name]['jacobian'] = component_jacobian[ix]/math.sqrt(2.)
    return result


def block_scores(residual, component):
    ans = {k: float(4*np.dot(residual[a:b], residual[a:b])) for k, (a, b) in SLICES.items()}
    ans.update(P_energy=float(np.mean(component[[0, 2]]**2)), B_energy=float(np.mean(component[[1, 3]]**2)))
    return ans


def comparisons(actual, origin):
    def nonincrease(k):
        return actual[k] <= origin[k]+1e-10*max(1., abs(actual[k]), abs(origin[k]))
    strict = {k: actual[k] < origin[k]-1e-10*max(1., abs(actual[k]), abs(origin[k])) for k in actual}
    return {'strict_descent_by_block': strict, 'nonincrease_by_block': {k: nonincrease(k) for k in actual},
            'both_component_scores_decrease': strict['P_energy'] and strict['B_energy'],
            'old_EF_nonincrease': nonincrease('old_EF'),
            'all_original_four_nonincrease': all(nonincrease(k) for k in SLICES)}


def load_context(inputs_path):
    inputs_pin = pin(inputs_path)
    src = json.loads(Path(inputs_path).read_text())
    data = {k: read(v) for k, v in src['artifacts'].items() if Path(v['path']).suffix == '.json'}
    for record in src['artifacts'].values():
        require(pin(record['path']) == record, 'Direct source pin')
    fit, review, root = (data[k] for k in ('fit_result', 'full_fit_review', 'root_fit_review'))
    require(fit['status'] == 'bounded_joint_response_fit_complete' and fit['error'] is None and fit['candidate_frozen'] is True, 'Complete frozen fit')
    require(len(fit['vectors']) == 673 and len(fit['iterations']) == 24 and fit['completed_outer_steps'] == 24, '673 vectors/24 origins')
    require(fit['vectors'][fit['selected_key']]['artifact'] == fit['selected'] == src['artifacts']['selected'], 'Exact selected artifact pin')
    require(review['status'] == 'passed_terminal_joint_fit_review_continuation' and review['aggregate_audit_status'] == 'passed_using_preserved_completed_phases_and_narrow_continuation', 'Actual continuation admission')
    require(root['status'] == 'admitted_completed_frozen_fit_provenance_and_numerical_review' and root['independent_review'] == src['artifacts']['full_fit_review'], 'Root full-fit admission')
    for record in (review, root):
        require(record['fit_result'] == src['artifacts']['fit_result'] and record['selected'] == fit['selected'], 'Admitted fit/selected join')
    em, dec, cfg = (data[k] for k in ('evaluator_manifest', 'decomposition_result', 'evaluator_configuration'))
    fc = data['fit_configuration']
    require(fit['configuration'] == src['artifacts']['fit_configuration'] and fit['implementation'] == src['artifacts']['fit_implementation'] and fc['evaluator_config'] == src['artifacts']['evaluator_configuration'] and fc['parameter_domain'] == src['artifacts']['domain'], 'Executed driver/config/evaluator lineage')
    require(cfg['manifest'] == src['artifacts']['evaluator_manifest'] and em['sources']['target_inventory'] == src['artifacts']['targets'], 'Evaluator static/target lineage')
    require(fit['parameter_domain'] == src['artifacts']['domain'] and em['sources']['domain'] == src['artifacts']['domain'], 'Exact domain lineage')
    require(dec['status'] == 'complete_saved_pair_collective_diagnostic' and dec['sources']['selected'] == fit['selected'], 'Selected decomposition lineage')
    for key in ('evaluator_manifest', 'targets', 'common_basis', 'matched_pair'):
        require(dec['sources'][key] == src['artifacts'][key], 'Decomposition direct source '+key)
    require(len(em['rows']) == 46 and len(em['static']) == 48, '46 variable+48 static rows')
    params = data['domain']['parameters']
    lower = np.asarray([p['lower'] for p in params]); upper = np.asarray([p['upper'] for p in params])
    require(len(params) == 12 and np.all(upper > lower), 'Twelve finite box variables')
    terms, refs = cycles(dec, data['targets'])
    return inputs_pin, src, data, lower, upper, terms, refs


def calculate(inputs_path, release_path, output):
    tick, cpu = time.monotonic(), time.process_time()
    require(all(os.environ.get(k) == '1' for k in THREADS), 'One numerical thread')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    result = {'status': 'started', 'implementation': pin(__file__), 'error': None,
              'new_native_calls': 0, 'new_reference_calls': 0, 'new_parameter_evaluations': 0,
              'new_fit_or_selection': False, 'fresh_outcomes_read': False}
    try:
        inputs_pin, src, data, lower, upper, terms, refs = load_context(inputs_path)
        require(pin(HERE/'LP.py') == src['lp'], 'Frozen LP source')
        require(pin(src['protocol']['path']) == src['protocol'], 'Frozen protocol')
        release = json.loads(Path(release_path).read_text())
        require(release['status'] == 'released_saved_joint_component_feasibility' and release['implementation'] == pin(__file__) and release['inputs'] == inputs_pin and release['protocol'] == src['protocol'] and release['lp'] == src['lp'] and release['native_program_cap'] == 0 and release['new_parameter_vector_cap'] == 0, 'Exact root saved-only release')
        spec = importlib.util.spec_from_file_location('frozen_component_lp', HERE/'LP.py')
        lp = importlib.util.module_from_spec(spec); spec.loader.exec_module(lp)
        fit, em, dec = (data[k] for k in ('fit_result', 'evaluator_manifest', 'decomposition_result'))
        static = {k: r['energy_kcal_mol'] for k, r in em['static'].items()}
        width = upper-lower
        vectors, residuals, components, parameters = {}, {}, {}, {}
        maximum_closure = 0.
        for key, record in fit['vectors'].items():
            require(record['status'] == 'complete', 'Every cached vector complete')
            v = read(record['artifact'])
            require(v['status'] == 'complete' and v['parameters'] == record['parameters'] and token(v['parameters']) == key, 'Exact vector/parameter identity')
            require(set(v['cells']) == set(em['rows']) and all(x['status'] == 'complete' for x in v['cells'].values()), 'Full46 cell denominator')
            z = (np.asarray(v['parameters'])-lower)/width
            require(np.all(np.isfinite(z)) and np.all(z >= -1e-12) and np.all(z <= 1+1e-12), 'Parameter domain')
            score = v['score']; r = np.asarray(score['weighted_residuals'], dtype=float)
            require(r.shape == (531,) and np.all(np.isfinite(r)), 'Finite original531 residuals')
            require({k: score['residual_layout'][v] for k, v in SCORE_KEYS.items()} == {k: list(v) for k, v in SLICES.items()}, 'Original residual ordering')
            energies = dict(static, **{k: r['energy'] for k, r in v['cells'].items()})
            require(all(math.isfinite(e) for e in energies.values()), 'Finite saved energies')
            parents, paired, e, closure = component_values(energies, terms, refs)
            maximum_closure = max(maximum_closure, closure)
            scores = block_scores(r, e)
            for name, old in SCORE_KEYS.items():
                close(scores[name], score['block_objectives'][old], 'Stored block score')
            close(float(r@r), score['objective'], 'Stored original J')
            allgates = gates(score)
            residuals[key], components[key], parameters[key] = r, e, np.asarray(v['parameters'])
            vectors[key] = {'artifact': record['artifact'], 'parameters': v['parameters'], 'role': record['role'],
                            'outer_step': record.get('outer_step'), 'parents_kcal_mol': parents, 'paired': paired,
                            'component_errors_kcal_mol': dict(zip(ERROR_ORDER, e.tolist())), 'scores': scores,
                            'original_objective': score['objective'], 'original25': allgates,
                            'original25_pass': score['original_EF_checks_passed']}
        for key, label in ((fit['initial_key'], 'baseline'), (fit['selected_key'], 'selected')):
            for parent in PARENTS:
                close(list(vectors[key]['parents_kcal_mol'][parent].values()), list(dec['parents'][parent][label+'_kcal_mol'].values()), 'Admitted '+label+' parent decomposition')
            for n in (8, 9):
                for channel in ('C', 'P', 'B'):
                    close(vectors[key]['paired']['W'+str(n)]['error_kcal_mol'][channel], dec['paired_Dy_minus_La']['W'+str(n)]['errors'][label]['model_minus_reference_kcal_mol'][channel], 'Admitted paired decomposition')
        origins, trials = [], []
        max_j_error = 0.
        for expected_step, entry in enumerate(fit['iterations'], 1):
            require(entry['outer_step'] == expected_step and entry['status'] == 'complete', 'Ordered complete24 origins')
            key = entry['origin_key']; center = parameters[key]; z = (center-lower)/width
            require(len(entry['derivative_keys']) == 24 and len(entry['trial_keys']) == len(entry['trial_points']) == 4, '24 probes/four trials per origin')
            cols, extra, spans = [], [], []
            for j in range(12):
                ka, kb = entry['derivative_keys'][2*j:2*j+2]
                for probe, sign in ((ka, -1), (kb, 1)):
                    target = z.copy(); target[j] = np.clip(target[j]+sign*.005, 0, 1)
                    expected = center.copy(); expected[j] = np.clip(center[j]+(target[j]-z[j])*width[j], lower[j], upper[j])
                    require(parameters[probe].tolist() == expected.tolist(), 'Exact probe parameter/order')
                span = (parameters[kb][j]-parameters[ka][j])/width[j]
                require(span > 0 and span == entry['derivative_spans'][j], 'Exact actual normalized probe span')
                spans.append(span); cols.append((residuals[kb]-residuals[ka])/span); extra.append((components[kb]-components[ka])/span)
            jac = np.column_stack(cols); jc = np.column_stack(extra)
            replay = close(jac, entry['Jacobian_normalized'], 'Saved531 Jacobian replay'); max_j_error = max(max_j_error, replay)
            rb = blocks(residuals[key], components[key], jac, jc)
            local = lp.solve(rb, z)
            origins.append({'outer_step': expected_step, 'key': key, 'original25_pass': vectors[key]['original25_pass'],
                            'normalized_parameters': z.tolist(), 'derivative_keys': entry['derivative_keys'], 'derivative_spans': spans,
                            'Jacobian_original531': jac.tolist(), 'Jacobian_components_order': list(ERROR_ORDER),
                            'Jacobian_components': jc.tolist(), 'original_Jacobian_max_abs_difference': replay,
                            'saved_linear_rank': entry['linear_directions'].get('rank'), 'LP': local})
            for index, (trial_key, point) in enumerate(zip(entry['trial_keys'], entry['trial_points'])):
                require(parameters[trial_key].tolist() == point['parameters'], 'Exact saved trial parameters')
                dz = (parameters[trial_key]-center)/width
                require(np.max(np.abs(dz)) <= .025+1e-12, 'Saved trial trust domain')
                predicted_r = residuals[key]+jac@dz; predicted_e = components[key]+jc@dz
                close(float(predicted_r@predicted_r), point['predicted_objective'], 'Saved trial predicted original J')
                actual_b = blocks(residuals[trial_key], components[trial_key]); pred_b = blocks(predicted_r, predicted_e)
                discrepancy = {k: {'residual_difference_L2': float(np.linalg.norm(actual_b[k]['residual']-pred_b[k]['residual'])),
                                   'residual_difference_max_abs': float(np.max(abs(actual_b[k]['residual']-pred_b[k]['residual']))),
                                   'actual_norm': float(np.linalg.norm(actual_b[k]['residual'])),
                                   'predicted_norm': float(np.linalg.norm(pred_b[k]['residual']))} for k in rb}
                trials.append({'outer_step': expected_step, 'trial_index': index, 'origin_key': key, 'trial_key': trial_key,
                               'method': point['method'], 'fraction': point['fraction'], 'normalized_delta': dz.tolist(),
                               'origin_original25_pass': vectors[key]['original25_pass'], 'trial_original25_pass': vectors[trial_key]['original25_pass'],
                               'actual_scores': vectors[trial_key]['scores'], 'predicted_scores': block_scores(predicted_r, predicted_e),
                               'actual_component_errors_kcal_mol': components[trial_key].tolist(), 'predicted_component_errors_kcal_mol': predicted_e.tolist(),
                               'actual_original_residuals531': residuals[trial_key].tolist(), 'predicted_original_residuals531': predicted_r.tolist(),
                               'actual_original25': vectors[trial_key]['original25'], 'linearization': discrepancy,
                               'changes': {k: vectors[trial_key]['scores'][k]-vectors[key]['scores'][k] for k in vectors[key]['scores']},
                               'flags': comparisons(vectors[trial_key]['scores'], vectors[key]['scores'])})
        require(len(vectors) == 673 and len(origins) == 24 and len(trials) == 96, 'Full diagnostic denominator')
        result.update(status='complete_saved_joint_component_feasibility', inputs=inputs_pin, sources=src,
                      root_release=pin(release_path), numerical_versions={'numpy': np.__version__, 'scipy': scipy.__version__},
                      component_order=list(ERROR_ORDER), reference_targets=refs, literal_cycles=terms,
                      vectors=vectors, origins=origins, trials=trials,
                      counts={'vectors': 673, 'origins': 24, 'saved_trials': 96, 'expected_nondegenerate_LPs': 72,
                              'LPs': sum(len(r['LP']['modes']) for r in origins)},
                      degenerate_origin_steps=[r['outer_step'] for r in origins if r['LP']['status'] == 'degenerate_block_norm'],
                      initial_key=fit['initial_key'], selected_key=fit['selected_key'], last_Jacobian_origin=fit['last_Jacobian_origin_key'],
                      last_origin_is_selected=fit['last_Jacobian_origin_key'] == fit['selected_key'],
                      origin_steps_passing_all25=[r['outer_step'] for r in origins if r['original25_pass']],
                      maximum_cycle_closure=maximum_closure, maximum_original_Jacobian_replay_error=max_j_error,
                      limits=['First-order aggregate LP does not preserve individual25 gates.', 'No selected-point Jacobian is inferred.',
                              'No global no-go, physical acceptance, new candidate selection or fresh validation.'])
    except Exception as exc:
        result.update(status='failed_saved_joint_component_feasibility', error=repr(exc))
    result['cost'] = {'wall_seconds': time.monotonic()-tick, 'CPU_seconds': time.process_time()-cpu,
                      'threads': 1, 'native_programs': 0, 'new_reference_calls': 0, 'scheduler_allocated_CPU_seconds': 0}
    save(out/'RESULT.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--inputs', required=True); p.add_argument('--root-release', required=True); p.add_argument('--out', required=True)
    a = p.parse_args(); r = calculate(a.inputs, a.root_release, a.out)
    print(json.dumps({k: r.get(k) for k in ('status', 'error', 'counts', 'cost')}, indent=2))
    raise SystemExit(0 if r['status'] == 'complete_saved_joint_component_feasibility' else 1)
