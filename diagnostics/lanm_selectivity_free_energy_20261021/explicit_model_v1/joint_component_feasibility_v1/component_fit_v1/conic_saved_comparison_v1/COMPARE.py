#!/usr/bin/env python3
"""One terminal-gated saved SOCP comparison. No evaluator or native execution."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time

import clarabel
import numpy as np
from scipy import sparse

THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
           'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')
BLOCKS = ('old_EF', 'pair_mu', 'C_mu', 'chi', 'P_energy', 'B_energy')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(p):
    require(isinstance(p, dict) and pin(p['path']) == p, 'Changed or unbound pin')
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    os.replace(temporary, path)


def module(p):
    spec = importlib.util.spec_from_file_location('frozen_component_core', verify(p))
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def attributes(value):
    return {k: getattr(value, k) for k in dir(value)
            if not k.startswith('_') and not callable(getattr(value, k))}


def plain(value):
    """Preserve nonfinite native diagnostics explicitly without invalid JSON."""
    if isinstance(value, dict):
        return {k: plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [plain(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else repr(float(value))
    if value is None or isinstance(value, (str, int, bool)):
        return value
    return str(value)


def fingerprint(config):
    return digest({k: v for k, v in config.items()
                   if k not in ('terminal_campaign', 'source_review', 'root_release')})


def frozen_settings(p):
    declared = read(p)
    settings = clarabel.DefaultSettings()
    require(set(attributes(settings)) == set(declared), 'All effective settings pinned')
    for key, value in declared.items():
        setattr(settings, key, value)
    require(attributes(settings) == declared, 'Exact effective settings')
    require(declared['max_iter'] == 200 and declared['time_limit'] == 60
            and declared['max_threads'] == 1 and declared['direct_solve_method'] == 'qdldl', 'Fixed solver resource limits')
    require(all(declared[k] == 1e-10 for k in
                ('tol_gap_abs', 'tol_gap_rel', 'tol_feas', 'tol_infeas_abs', 'tol_infeas_rel')), 'Fixed primary tolerances')
    return settings


def terminal_gate(fit, cfg, inputs):
    require(fit['status'] in ('bounded_component_fit_complete', 'failed_component_fit')
            and fit['mode'] == 'run' and 'cost' in fit, 'Campaign must already be terminal')
    require(fit['implementation'] == inputs['campaign_driver']
            and fit['configuration'] == inputs['campaign_configuration']
            and fit['configuration_fingerprint'] == inputs['campaign_fingerprint'], 'Exact current campaign only')
    require(cfg['terminal_campaign']['path'] == inputs['expected_terminal_path'], 'Literal terminal result path')


def context(cfg):
    require(cfg['implementation'] == pin(__file__), 'Exact comparison source')
    require(all(os.environ.get(k) == '1' for k in THREADS), 'One numerical thread')
    inputs = read(cfg['inputs'])
    env = read(cfg['environment'])
    require(str(Path(sys.prefix).resolve()) == env['prefix'], 'Isolated environment only')
    require(clarabel.__version__ == '0.11.1', 'Declared Clarabel release')
    import importlib.metadata as metadata
    for name, version in env['versions'].items():
        require(metadata.version(name) == version, 'Pinned package version ' + name)
    for p in env['files']:
        verify(p)
    verify(env['installation']); verify(env['pip_report'])
    settings = frozen_settings(cfg['settings'])
    # No scientific vectors or partial campaign outcomes are read before this gate.
    fit = read(cfg['terminal_campaign'])
    terminal_gate(fit, cfg, inputs)
    review = read(cfg['source_review']); release = read(cfg['root_release'])
    require(review['status'] == 'passed_conic_saved_comparison_source_review'
            and review['implementation'] == cfg['implementation']
            and review['configuration_fingerprint'] == fingerprint(cfg), 'Exact source review')
    require(release['status'] == 'released_one_saved_origin24_conic_comparison'
            and release['configuration_fingerprint'] == fingerprint(cfg)
            and release['terminal_campaign'] == cfg['terminal_campaign']
            and release['source_review'] == cfg['source_review']
            and release['maximum_solver_calls'] == 1 and release['new_native_calls'] == 0, 'One root-released comparison')
    return inputs, fit, settings


def load_problem(inputs, core):
    campaign_cfg = read(inputs['campaign_configuration'])
    require(campaign_cfg['core'] == inputs['core'] and campaign_cfg['reuse_inventory'] == inputs['reuse_inventory'], 'Frozen campaign source/imports')
    inv = read(inputs['reuse_inventory']); declaration = read(inputs['campaign_declaration'])
    root = read(inv['root_fit_review']); audit = read(root['independent_review'])
    require(root['status'] == 'admitted_completed_frozen_fit_provenance_and_numerical_review'
            and audit['status'] == 'passed_terminal_joint_fit_review_continuation'
            and root['fit_result'] == inv['fit_result'] == audit['fit_result'], 'Original saved fit admission')
    verify(inv['fit_result'])
    require(inv['origin_key'] == declaration['start']['key'] == inputs['origin_key']
            and inv['vectors'][inv['origin_key']] == declaration['start']['artifact'], 'Exact declared origin24')
    require(inputs['derivative_keys'] == inv['derivative_keys'] and len(inputs['derivative_keys']) == 24
            and inputs['vectors'] == inv['vectors']
            and set(inputs['vectors']) == {inputs['origin_key'], *inputs['derivative_keys']}, 'Only immutable origin and first24 probes')
    values = {k: read(p) for k, p in inputs['vectors'].items()}
    for key, value in values.items():
        require(value['status'] == 'complete' and core.token(value['parameters']) == key, 'Exact admitted vector identity')
    manifest = read(inputs['manifest']); domain = read(inputs['domain'])
    parent = read(campaign_cfg['parent_configuration'])
    require(parent['manifest'] == inputs['manifest'] and parent['domain'] == inputs['domain']
            and campaign_cfg['decomposition'] == inputs['decomposition']
            and campaign_cfg['declaration'] == inputs['campaign_declaration'], 'Same scientific authorities')
    mapping_review = read(inputs['mapping_review'])
    require(mapping_review['status'] == 'passed_exact_conic_formulation_source_review'
            and mapping_review['implementation'] == inputs['core'] and mapping_review['mapping'] == inputs['mapping'], 'Reviewed mapping source')
    verify(inputs['mapping']); verify(inputs['declaration'])
    origin = values[inputs['origin_key']]
    fixed = core.make_fixed(origin, manifest, read(inputs['decomposition']), domain, declaration)
    problem = core.linearize(origin, [values[k] for k in inputs['derivative_keys']], fixed)
    return problem, fixed


def assemble(p):
    """Clarabel form A*x+s=b, s in K, x=(d_0..d_11,t)."""
    require(tuple(p['blocks']) == BLOCKS and len(p['forces']) == 15 and len(p['energy']['r']) == 10, 'Original six blocks and25 gates')
    rows, rhs, cones, layout = [], [], [], []

    def soc(name, r, J, scale, component=False):
        start = len(rhs)
        top = np.zeros(13)
        if component:
            top[-1] = -1
        rows.append(top); rhs.append(0. if component else 1.)
        for value, row in zip(np.asarray(r) * scale, np.asarray(J) * np.asarray(scale)[..., None] if np.ndim(scale) else np.asarray(J) * scale):
            rows.append(np.r_[-row, 0.]); rhs.append(float(value))
        size = len(rhs) - start
        cones.append(clarabel.SecondOrderConeT(size))
        layout.append({'name': name, 'kind': 'SOC', 'start': start, 'stop': len(rhs)})

    for name, block in p['blocks'].items():
        require(block['denominator'] > 0, 'Positive normalization')
        soc('block/' + name, block['r'], block['J'], 1/math.sqrt(block['denominator']), name in ('P_energy', 'B_energy'))
    for name, block in p['forces'].items():
        require(np.all(block['weights'] > 0), 'Positive original RMS weights')
        soc('force/RMS/' + name, block['r'], block['J'], np.sqrt(block['weights'])/2)

    linear_start = len(rhs)

    def side(name, constant, gradient):
        rows.append(-np.asarray(gradient, dtype=float)); rhs.append(float(constant))
        layout.append({'name': name, 'kind': 'nonnegative', 'start': len(rhs)-1, 'stop': len(rhs)})

    for name, r, J in zip(p['energy']['names'], p['energy']['r'], p['energy']['J']):
        side('energy/upper/' + name, 1-r, np.r_[-J, 0.])
        side('energy/lower/' + name, 1+r, np.r_[J, 0.])
    for name, block in p['forces'].items():
        for i, (r, J) in enumerate(zip(block['r'], block['J'])):
            side('force/upper/' + name + '/' + str(i), 1-r/5, np.r_[-J/5, 0.])
            side('force/lower/' + name + '/' + str(i), 1+r/5, np.r_[J/5, 0.])
    lower = np.r_[p['bounds_lower'][:12], 0.]
    upper = np.r_[p['bounds_upper'][:12], math.sqrt(p['q_current'])]
    for i in range(13):
        row = np.eye(13)[i]
        side('bound/lower/' + str(i), -lower[i], row)
        side('bound/upper/' + str(i), upper[i], -row)
    cones.append(clarabel.NonnegativeConeT(len(rhs)-linear_start))
    A = sparse.csc_matrix(np.asarray(rows))
    b = np.asarray(rhs); c = np.r_[np.zeros(12), 1.]
    require(linear_start == 961 and len(rhs)-linear_start == 856 and len(cones) == 22 and A.shape == (1817, 13), 'Exact21 SOC and856 linear sides')
    require(np.all(np.isfinite(A.data)) and np.all(np.isfinite(b)), 'Finite conic problem')
    return {'P': sparse.csc_matrix((13, 13)), 'c': c, 'A': A, 'b': b, 'cones': cones, 'layout': layout}


def correspondence(p, model, x, core):
    original, _, names = core.constraints(p, np.r_[x[:12], x[-1]**2])
    slack = model['b'] - model['A'] @ x
    reconstructed = {}
    for row in model['layout']:
        s = slack[row['start']:row['stop']]
        if row['kind'] == 'SOC':
            reconstructed[row['name']] = float(s[0]**2 - np.dot(s[1:], s[1:]))
        elif not row['name'].startswith('bound/'):
            reconstructed[row['name']] = float(s[0])
    require(set(reconstructed) == set(names) and len(names) == 851, 'Complete original constraint names')
    replay = np.asarray([reconstructed[n] for n in names])
    # Arithmetic identity only; never substitutes for exact physical certificates.
    core.near(replay, original, 'Cone/CORE affine and squared-norm correspondence')
    return {'original_constraint_count': len(names), 'maximum_arithmetic_difference': float(np.max(abs(replay-original))),
            'nonnegative_t': bool(x[-1] >= 0), 'minimum_core_margin': float(np.min(original))}


def cone_margins(model, vector):
    result = {}
    for row in model['layout']:
        s = vector[row['start']:row['stop']]
        result[row['name']] = float(s[0]-np.linalg.norm(s[1:])) if row['kind'] == 'SOC' else float(s[0])
    return result


def diagnostics(model, solution):
    x, s, z = (np.asarray(getattr(solution, k), dtype=float) for k in ('x', 's', 'z'))
    if not all(np.all(np.isfinite(v)) for v in (x, s, z)):
        return {'finite_return': False, 'certificate': 'unavailable_nonfinite_return'}
    primal = model['A'] @ x + s - model['b']
    dual = model['A'].T @ z + model['c']
    pcost = float(model['c'] @ x); dcost = float(-model['b'] @ z)
    return {'finite_return': True, 'primal_equation_residual': primal.tolist(), 'dual_stationarity_residual': dual.tolist(),
            'primal_residual_inf': float(np.max(abs(primal))), 'dual_residual_inf': float(np.max(abs(dual))),
            'primal_objective_t': pcost, 'dual_objective_t': dcost, 'signed_gap_t': pcost-dcost,
            'complementarity_dot': float(s @ z), 'primal_affine_cone_margins': cone_margins(model, model['b']-model['A']@x),
            'raw_slack_cone_margins': cone_margins(model, s), 'dual_cone_margins': cone_margins(model, z),
            'rigorously_rounded_dual_certificate': False}


def fraction_reports(p, x, status, core):
    require(len(x) == 13 and np.all(np.isfinite(x)) and x[-1] >= 0, 'Finite return with nonnegative t before squaring')
    returned = np.r_[x[:12], x[-1]**2]
    slots = [core.certify_fraction(p, returned, f, status == 'Solved') for f in core.FRACTIONS]
    epsilon = 1e-10 * max(1., abs(p['q_current']))
    for slot in slots:
        slot['q_selection_epsilon'] = epsilon
        slot['predicted_q_improvement_exceeds_selection_epsilon'] = slot['certificate']['q_fraction'] < p['q_current']-epsilon
        slot['native_evaluated'] = False
    return slots


def run(config_path):
    tick, cpu = time.monotonic(), time.process_time()
    cfg = json.loads(Path(config_path).read_text())
    inputs, fit, settings = context(cfg)
    out = Path(cfg['output_directory']); out.mkdir(parents=True, exist_ok=False)
    result = {'status': 'started', 'implementation': pin(__file__), 'configuration': pin(config_path),
              'configuration_fingerprint': fingerprint(cfg), 'terminal_campaign': cfg['terminal_campaign'],
              'solver_calls': 0, 'new_native_calls': 0, 'new_candidate_evaluations': 0, 'error': None,
              'selected_candidate': None, 'automatic_fit_extension': False, 'original_SLSQP_modified': False}
    save(out/'RESULT.json', result)
    try:
        core = module(inputs['core']); p, fixed = load_problem(inputs, core); model = assemble(p)
        result.update(origin_key=inputs['origin_key'], derivative_keys=inputs['derivative_keys'],
                      fixed_ceilings=fixed['ceilings'], fixed_denominators=fixed['denominators'],
                      derivative_spans=p['spans'], q_current=p['q_current'], effective_settings=attributes(settings))
        result['zero_step_correspondence'] = correspondence(p, model, np.r_[np.zeros(12), math.sqrt(p['q_current'])], core)
        zero, _, _ = core.constraints(p, np.r_[np.zeros(12), p['q_current']])
        require(float(np.min(zero)) >= -1e-10, 'Inherited CORE zero-step startup arithmetic gate')
        # Read only the original first SLSQP proposal after terminal status. Never later origins/candidates.
        if fit['iterations'] and fit['iterations'][0].get('proposal'):
            first = fit['iterations'][0]
            require(first['origin_key'] == inputs['origin_key'] and first['derivative_keys'] == inputs['derivative_keys'], 'Original first proposal identity')
            prior = read(first['proposal'])
            require(prior['current_key'] == inputs['origin_key'] and prior['derivative_keys'] == inputs['derivative_keys'], 'Original SLSQP model identity')
            core.near(prior['Jacobian_original531'], p['Jacobian_original531'], 'Original531 Jacobian')
            core.near(prior['Jacobian_components'], p['Jacobian_components'], 'Original component Jacobian')
            core.near(prior['derivative_spans'], p['spans'], 'Original actual normalized probe spans')
            require(prior['fixed_ceilings'] == fixed['ceilings'] and prior['fixed_component_denominators'] == fixed['denominators'], 'Original normalizations')
            result['original_SLSQP'] = {'artifact': first['proposal'], 'status': prior['status'], 'solver': prior.get('solver'),
                                       'fraction_eligibility': [{'fraction': s['fraction'], 'eligible': s['eligible']} for s in prior['trial_slots']]}
        else:
            result['original_SLSQP'] = {'status': 'no_first_proposal_in_terminal_campaign'}
        np.savez(out/'CONIC_DATA.npz', A_data=model['A'].data, A_indices=model['A'].indices, A_indptr=model['A'].indptr,
                 A_shape=model['A'].shape, b=model['b'], c=model['c'])
        save(out/'ROW_LAYOUT.json', model['layout'])
        save(out/'LINEARIZED_PROBLEM.json', plain(p))
        result['linearized_problem'] = pin(out/'LINEARIZED_PROBLEM.json')
        result['conic_data'] = pin(out/'CONIC_DATA.npz'); result['row_layout'] = pin(out/'ROW_LAYOUT.json')
        # A fixed exclusive marker survives crashes and forbids duplicate attempts across output versions.
        with open(cfg['one_attempt_marker'], 'x') as stream:
            json.dump({'configuration': pin(config_path), 'terminal_campaign': cfg['terminal_campaign'], 'pid': os.getpid()}, stream)
        result['attempt_reserved'] = True; save(out/'RESULT.json', result)
        solver = clarabel.DefaultSolver(model['P'], model['c'], model['A'], model['b'], model['cones'], settings)
        require(attributes(solver.get_settings()) == attributes(settings), 'Backend effective settings unchanged')
        solver.print_to_buffer()
        result['solver_calls'] = 1; save(out/'RESULT.json', result)
        solution = solver.solve()  # Sole solver call. No retry or fallback.
        info = solver.get_info()
        (out/'SOLVER.log').write_text(solver.get_print_buffer())
        result['solution'] = plain(attributes(solution))
        raw_info = attributes(info); raw_info['linsolver'] = attributes(info.linsolver)
        result['info'] = plain(raw_info)
        result['solver_status'] = str(solution.status)
        result['solver_success'] = result['solver_status'] == 'Solved'
        result['reduced_accuracy_status'] = result['solver_status'].startswith('Almost')
        # Persist raw return before any optional independent arithmetic can fail.
        save(out/'RESULT.json', result)
        result['independent_linear_residuals'] = plain(diagnostics(model, solution))
        x = np.asarray(solution.x)
        if x.shape == (13,) and np.all(np.isfinite(x)) and x[-1] >= 0:
            result['returned_q'] = float(x[-1]**2)
            result['returned_correspondence'] = correspondence(p, model, x, core)
            result['fractions'] = fraction_reports(p, x, result['solver_status'], core)
        else:
            result['fractions'] = []; result['fraction_rejection'] = 'Nonfinite, malformed or negative-t return; no clipping/squaring repair'
        result['status'] = 'completed_saved_conic_comparison'
        result['interpretation'] = 'Solver status, floating cone residuals, and exact source-residual fraction certificates are separate. No native candidate, physical validation, or globally optimal nonlinear model is established.'
    except Exception as exc:
        result.update(status='failed_saved_conic_comparison', error=repr(exc))
    finally:
        result['cost'] = {'wall_seconds': time.monotonic()-tick, 'CPU_seconds': time.process_time()-cpu, 'numerical_threads': 1, 'scheduler_job_id': None}
        save(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json'), 'solver_calls': result['solver_calls'], 'error': result['error']}))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--config', required=True, type=Path)
    arguments = parser.parse_args()
    answer = run(arguments.config)
    raise SystemExit(bool(answer['error']))
