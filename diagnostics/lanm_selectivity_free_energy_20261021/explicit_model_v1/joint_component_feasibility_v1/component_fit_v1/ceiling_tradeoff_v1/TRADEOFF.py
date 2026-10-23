#!/usr/bin/env python3
"""One separately released saved-matrix tradeoff; no native/candidate API."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time

import clarabel
import numpy as np

import BASELINE as baseline

require, pin, read = baseline.require, baseline.pin, baseline.read


def fingerprint(cfg):
    v = {k: x for k, x in cfg.items() if k not in ('source_review', 'root_release')}
    return hashlib.sha256(json.dumps(v, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def witness(p, returned, fraction, solver_success, core, removed, reference_rho):
    slot = core.certify_fraction(p, returned, fraction, solver_success)
    raw = slot['certificate']
    rho, cushion, q = baseline.auxiliary_q([raw['blocks'][k]['ratio'] for k in ('P_energy', 'B_energy')])
    # Re-express the exact physical fraction as fraction one, replacing only q.
    # This avoids applying fraction interpolation to the new auxiliary headroom.
    constructed = core.certify_fraction(p, np.r_[fraction*np.asarray(returned[:12]), q], 1., solver_success)
    require(constructed['parameters'] == slot['parameters'] and constructed['key'] == slot['key'], 'Witness physical fraction unchanged')
    c = constructed['certificate']
    retained = {k: v for k, v in c['residual_constraints'].items() if k not in removed}
    require(len(retained) == 847, 'Exactly847 retained scalar constraints')
    feasible = bool(c['in_box_and_trust'] and c['all25_predicted']
                    and all(math.isfinite(v) and v > 0 for v in retained.values()))
    require(rho == max(c['blocks'][k]['ratio'] for k in ('P_energy', 'B_energy')), 'Unchanged actual component ratios')
    epsilon = 1e-10 * max(1., abs(reference_rho))
    return {'fraction': fraction, 'raw_source_certificate': slot, 'constructed_source_certificate': constructed,
            'raw_auxiliary_q': raw['q_fraction'], 'raw_minimum_epigraph_margin': raw['minimum_epigraph_margin'],
            'auxiliary_cushion': cushion, 'physical_fraction_changed': False,
            'retained_scalar_margin_count': len(retained), 'minimum_retained_margin': min(retained.values()),
            'strict_retained_feasible': feasible, 'actual_component_rho': rho,
            'auxiliary_q': c['q_fraction'], 'baseline_actual_component_rho': reference_rho,
            'gain_vs_baseline_actual_rho': reference_rho-rho, 'gain_comparison_epsilon': epsilon,
            'feasible_gain_exceeds_numerical_comparison': feasible and rho < reference_rho-epsilon,
            'removed_ceiling_ratios': {k: c['blocks'][k]['ratio'] for k in ('old_EF', 'pair_mu', 'C_mu', 'chi')},
            'native_evaluated': False, 'candidate_admitted': False, 'nonlinear_or_global_optimality_claim': False}


def context(cfg):
    require(cfg['implementation'] == pin(__file__), 'Exact new diagnostic source')
    inputs = read(cfg['inputs'])
    require(pin(cfg['witness_contract']['path']) == cfg['witness_contract'], 'Prospective auxiliary-witness amendment')
    require(inputs['baseline_implementation'] == pin(baseline.__file__), 'Exact baseline helper')
    require(cfg['output_directory'] == inputs['solve_output_directory']
            and cfg['one_attempt_marker'] == inputs['one_attempt_marker'], 'Fixed output and one-attempt ledger')
    own = read(cfg['source_review']); release = read(cfg['root_release'])
    require(own['status'] == 'passed_ceiling_tradeoff_source_review'
            and own['implementation'] == cfg['implementation'] and own['inputs'] == cfg['inputs']
            and own['configuration_fingerprint'] == fingerprint(cfg), 'Exact source admission')
    require(release['status'] == 'released_one_saved_ceiling_tradeoff'
            and release['configuration_fingerprint'] == fingerprint(cfg)
            and release['source_review'] == cfg['source_review']
            and release['baseline_certificate'] == cfg['baseline_certificate']
            and release['maximum_solver_calls'] == 1 and release['new_native_calls'] == 0, 'One released attempt')
    cert = read(cfg['baseline_certificate'])
    require(cert['status'] == 'passed_saved_auxiliary_epigraph_baseline' and cert['error'] is None
            and cert['implementation'] == inputs['baseline_implementation'] and cert['inputs'] == cfg['inputs']
            and cert['solver_calls'] == 0 and cert['native_calls'] == 0, 'Actual unchanged-vector certificate')
    base, old, core, p, model, full, settings = baseline.load(inputs)
    replay = baseline.certify(p, full, core)
    require(replay == cert['baseline'], 'Exact scalar baseline replay')
    return inputs, base, old, core, p, model, full, settings, cert


def run(config_path):
    wall, cpu = time.monotonic(), time.process_time()
    cfg = json.loads(Path(config_path).read_text())
    # Admission failures happen before an attempt is reserved or a solver exists.
    inputs, base, old, core, p, full_model, full, settings, cert = context(cfg)
    out = Path(cfg['output_directory']); out.mkdir(parents=True, exist_ok=False)
    result = {'status': 'started', 'implementation': pin(__file__), 'configuration': pin(config_path),
              'configuration_fingerprint': fingerprint(cfg), 'inputs': cfg['inputs'],
              'baseline_certificate': cfg['baseline_certificate'], 'full_comparison': inputs['full_result'],
              'solver_calls': 0, 'new_native_calls': 0, 'new_candidate_evaluations': 0,
              'selected_candidate': None, 'original_gate_changed': False, 'error': None}
    base.save(out/'RESULT.json', result)
    try:
        model, removed, indices = old.remove_joint_ceilings(full_model)
        result.update(origin_key=full['origin_key'], derivative_keys=full['derivative_keys'],
                      linearized_problem=full['linearized_problem'], removed_cones=removed,
                      retained_original_row_indices=indices, effective_settings=base.attributes(settings),
                      baseline_actual_component_rho=cert['baseline']['actual_component_rho'],
                      baseline_auxiliary_q=cert['baseline']['auxiliary_q'])
        np.savez(out/'CONIC_DATA.npz', A_data=model['A'].data, A_indices=model['A'].indices,
                 A_indptr=model['A'].indptr, A_shape=model['A'].shape, b=model['b'], c=model['c'])
        base.save(out/'ROW_LAYOUT.json', model['layout'])
        result['conic_data'] = pin(out/'CONIC_DATA.npz'); result['row_layout'] = pin(out/'ROW_LAYOUT.json')
        # Exclusive persistent reservation survives failure; no alternate-output retry.
        with open(cfg['one_attempt_marker'], 'x') as stream:
            json.dump({'configuration': pin(config_path), 'baseline_certificate': cfg['baseline_certificate'],
                       'pid': os.getpid()}, stream)
            stream.flush(); os.fsync(stream.fileno())
        result['attempt_reserved'] = True; base.save(out/'RESULT.json', result)
        solver = clarabel.DefaultSolver(model['P'], model['c'], model['A'], model['b'], model['cones'], settings)
        require(base.attributes(solver.get_settings()) == base.attributes(settings), 'All effective settings unchanged')
        solver.print_to_buffer(); result['solver_calls'] = 1; base.save(out/'RESULT.json', result)
        solution = solver.solve()  # Exactly one call. No retry, fallback or new parameter evaluation.
        info = solver.get_info(); (out/'SOLVER.log').write_text(solver.get_print_buffer())
        result['solution'] = base.plain(base.attributes(solution))
        info_values = base.attributes(info); info_values['linsolver'] = base.attributes(info.linsolver)
        result['info'] = base.plain(info_values); result['solver_status'] = str(solution.status)
        result['solver_success'] = result['solver_status'] == 'Solved'
        base.save(out/'RESULT.json', result)  # Raw return survives later reporting errors.
        result['independent_linear_residuals'] = base.plain(base.diagnostics(model, solution))
        x = np.asarray(solution.x, dtype=float)
        if x.shape == (13,) and np.all(np.isfinite(x)) and x[-1] >= 0:
            result['returned_q'] = float(x[-1]**2)
            result['retained_correspondence'] = old.correspondence(p, model, x, core)
            result['witnesses'] = [witness(p, np.r_[x[:12], x[-1]**2], f, result['solver_success'],
                                          core, old.REMOVED, cert['baseline']['actual_component_rho'])
                                   for f in core.FRACTIONS]
        else:
            result['witnesses'] = []; result['witness_rejection'] = 'Malformed/nonfinite/negative epigraph; no repair'
        result['status'] = 'completed_saved_ceiling_tradeoff'
        result['interpretation'] = ('Status and numerical witness feasibility are separate. Feasible gain is a local '
            'tradeoff diagnostic; absent gain is inconclusive without adequate optimality evidence. No model/native admission.')
    except Exception as exc:
        result.update(status='failed_saved_ceiling_tradeoff', error=repr(exc))
    finally:
        result['cost'] = {'wall_seconds': time.monotonic()-wall, 'CPU_seconds': time.process_time()-cpu,
                          'numerical_threads': 1, 'scheduler_job_id': None}
        base.save(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json'),
                      'solver_calls': result['solver_calls'], 'error': result['error']}))
    return result


if __name__ == '__main__':
    a = argparse.ArgumentParser(); a.add_argument('--config', required=True)
    args = a.parse_args(); raise SystemExit(bool(run(args.config)['error']))
