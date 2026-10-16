#!/usr/bin/env python3
"""Matched 7/10-coordinate own-metal continuation; all numerical work on Slurm.

--stage preflight|scout|remaining --config CONFIG --output FRESH
Config pins implementation, plan, range_scout_result, range_remaining_result;
remaining additionally pins basis_scout_result. Physical POINT recipes stay
unchanged. All eight arm/state fixtures pass zero-call checks before model load.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import socket
import sys
import time
import traceback

PROTOCOL = 'lbt3_common_basis_response_v1'
RANGE_SHA = 'c61a12b1f4d47b984e0121c407c91f92e11288ae358e31f985969701b89c692b'
STATES = (('7CCO', 'La'), ('7CCO', 'Dy'), ('7CCN', 'La'), ('7CCN', 'Dy'))
ARMS = ('control', 'candidate')
EXTRA = (9, 38, 32)
EXTRA_NAMES = ('A3:sidechain:CB-CG', 'A11:sidechain:CB-CG', 'A10:phi:N-CA')
FD_STEPS = (1e-5, 5e-6)
FALLBACK = {'entry': 'all four initial trials rejected by geometry before any energy/cache evaluation',
    'coordinate_order': 'ascending numeric mode index', 'ranking': 'decreasing predicted own-energy decrease',
    'maximum_ranked_trials': 2, 'second_after_first_fallback_Armijo_rejection': True,
    'backend_failure': 'stop search unavailable', 'recursive_fallback': False, 'accepted_fallback_resets_history': True}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(record):
    path = Path(record['path']).resolve(strict=True)
    require(path.is_file() and pin(path)['sha256'] == record['sha256'], 'Missing/changed input: ' + record['path'])
    return path


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False); handle.write('\n')


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def original_metric(jac0, active, inherited=None):
    n = len(active); matrix = jac0[active].reshape(n, -1).T
    scales = np.linalg.norm(matrix, axis=0); require(np.all(scales > 1e-12), 'Zero source mode')
    normed = matrix/scales; eig, vectors = np.linalg.eigh(normed.T @ normed)
    require(int(np.sum(eig > 1e-10*eig[-1])) == n, 'Source arm metric rank deficient')
    transform = ((vectors/np.sqrt(eig)) @ vectors.T)/scales[:, None] if inherited is None else inherited.copy()
    require(transform.shape == (n, n) and np.allclose((matrix@transform).T@(matrix@transform), np.eye(n),
            rtol=0, atol=1e-10), 'Original-source whitening differs')
    return transform, {'source_column_norms': scales.tolist(), 'source_normalized_Gram_eigenvalues': eig.tolist(),
        'coordinate_scaling_T': transform.tolist(), 'source_metric_rank': n,
        'source_metric_condition': float(eig[-1]/eig[0]), 'metric_relative_cutoff': 1e-10}


def jacobian_check(kin, q):
    _, _, jac, _ = kin.evaluate(q); checks = []
    for k in EXTRA:
        for h in FD_STEPS:
            d = np.zeros(64); d[k] = h
            plus = kin.positions_only(q+d); minus = kin.positions_only(q-d)
            error = float(np.max(np.abs((plus-minus)/(2*h)-jac[k])))
            checks.append({'mode_index': k, 'mode_name': kin.modes[k]['name'], 'step_radian': h,
                'maximum_absolute_derivative_error_A_per_radian': error, 'pass': error <= 1e-6})
    require(all(r['pass'] for r in checks), 'Actual-fixture two-step Jacobian check failed')
    return checks


def inputs(config_path):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    cfg = json.loads(config_path.read_text()); verify(cfg['plan']); verify(cfg['implementation'])
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Common-basis adapter differs')
    scout, remaining = read(cfg['range_scout_result']), read(cfg['range_remaining_result'])
    require(scout['status'] == remaining['status'] == 'complete' and scout['stage'] == 'scout'
            and remaining['stage'] == 'remaining' and remaining['range_scout_result'] == cfg['range_scout_result']
            and scout['implementation']['sha256'] == remaining['implementation']['sha256'] == RANGE_SHA,
            'Complete unchanged range parents required')
    for parent in (scout, remaining):
        parent_cfg = read(parent['config']); verify(parent['implementation'])
        require(parent['implementation'] == parent_cfg['implementation'] and parent['plan'] == parent_cfg['plan']
                and parent['old_pool_result'] == parent_cfg['old_pool_result'] == scout['old_pool_result'], 'Range parent config differs')
    range_module = module('lbt3_basis_frozen_range', scout['implementation'])
    (range_cfg, old_pool, old_cfg, pool, base, world, old_active, settings, old_starts,
     old_checks, lbfgs_direction) = range_module.inputs(verify(scout['config']))
    global np
    import numpy as np
    require(scout['settings'] == remaining['settings'] == settings
            and scout['profile_result'] == remaining['profile_result'] == old_cfg['profile_result']
            and scout['fd_result'] == remaining['fd_result'] == old_cfg['fd_result'], 'Range recipe differs')
    from scaffold_accommodation import Geometry, SETTINGS as GEOMETRY_SETTINGS
    slots = [dict(row, parent=cfg[key]) for key, p in (('range_scout_result', scout), ('range_remaining_result', remaining)) for row in p['rows']]
    require(len(slots) == 4 and {(r['source_id'], r['metal']) for r in slots} == set(STATES), 'Range state denominator differs')
    active_by = {'control': old_active, 'candidate': old_active+list(EXTRA)}
    starts, preflight, metrics, fd_original, shared_cache, shared_records = {}, [], {}, {}, {}, {}
    for slot in sorted(slots, key=lambda r: STATES.index((r['source_id'], r['metal']))):
        source_id, metal = slot['source_id'], slot['metal']; state = (source_id, metal)
        require(read(slot['result']) == {k: v for k, v in slot.items() if k not in ('result', 'parent')}, 'Range row pin differs')
        old_start = old_starts[state]; source = old_start['source']; kin = source['kin']
        require(tuple(kin.modes[k]['name'] for k in EXTRA) == EXTRA_NAMES, 'Declared additional modes differ')
        search = read(slot['search'])
        require(slot['status'] == 'complete' and slot['candidate_admitted'] and search['candidate_admitted']
                and search['status'] == 'finite_proposal_available' and search['source_id'] == source_id and search['metal'] == metal
                and search['continuation_start'] == old_start['point_pin'] and search['previous_search'] == old_start['search_pin']
                and search['source_rebased'] is False and search['active_indices'] == old_active
                and search['active_names'] == base.NAMES and search['candidate'] == search['trace'][-1]['point']
                and search['candidate'] in slot['points'] and search['coordinate_scaling_T'] == old_start['search']['coordinate_scaling_T']
                and search['authoritative_origin'] == old_start['search']['authoritative_origin'], 'Range terminal lineage differs')
        cache, records = {}, []
        for record in slot['points']:
            point, _, _ = range_module.saved_point(pool, world, source, metal, record)
            if point['cache_key'] not in cache:
                cache[point['cache_key']] = record; records.append(record)
            else:
                require(read(cache[point['cache_key']])['components'] == point['components'], 'Conflicting prior point cache')
        point, x, gradient = range_module.saved_point(pool, world, source, metal, search['candidate'])
        q = np.asarray(search['full_q'], dtype=float); old_u = np.asarray(search['u_A'], dtype=float)
        require(q.shape == (64,) and np.isfinite(q).all() and np.count_nonzero(q[list(EXTRA)]) == 0
                and np.count_nonzero(q[[k for k in range(64) if k not in old_active]]) == 0
                and np.allclose(q[old_active], old_start['transform']@old_u, rtol=0, atol=1e-12), 'Range absolute q/old metric differs')
        replay, _, jac, _ = kin.evaluate(q); _, _, jac0, _ = kin.evaluate(np.zeros(64))
        require(np.allclose(replay, x, rtol=0, atol=1e-12), 'Range start kinematic replay differs')
        if source_id not in fd_original:
            fd_original[source_id] = jacobian_check(kin, np.zeros(64))
        fd_start = jacobian_check(kin, q)
        shared_cache[state], shared_records[state] = cache, records
        for arm in ARMS:
            active = active_by[arm]; metric_key = (source_id, arm)
            if metric_key not in metrics:
                metrics[metric_key] = original_metric(jac0, active, old_start['transform'] if arm == 'control' else None)
            transform, metric = metrics[metric_key]
            guard = range_module.range_guard(base, source, active, Geometry, GEOMETRY_SETTINGS, metal)
            origin_guard = guard.check(np.zeros(64), source['x'], source['x']); start_guard = guard.check(q, x, x)
            require(origin_guard['pass'] and start_guard['pass'], 'Arm original/start physical check failed')
            starts[(source_id, metal, arm)] = dict(source=source, search=search, search_pin=slot['search'],
                point=point, point_pin=search['candidate'], origin=old_start['origin'], origin_pin=old_start['origin_pin'],
                x=x, gradient=gradient, q=q, jac=jac, transform=transform, metric=metric, active=active,
                guard=guard, origin_guard=origin_guard, start_guard=start_guard, cache=cache, records=records)
            preflight.append({'source_id': source_id, 'metal': metal, 'arm': arm, 'pass': True,
                'previous_search': slot['search'], 'continuation_start': search['candidate'], 'original_point': old_start['origin_pin'],
                'active_indices': active, 'active_names': [kin.modes[k]['name'] for k in active], **metric,
                'initial_offset_v_A': [0.]*len(active), 'initial_geometry_unchanged': True,
                'start_replay_error_A': float(np.max(np.abs(replay-x))),
                'Jacobian_original': fd_original[source_id], 'Jacobian_range_start': fd_start,
                'Jacobian_probe_policy': 'centered mathematical map probes; no energy or claim that both signed probes lie inside search bounds',
                'original_guard': origin_guard, 'start_guard': start_guard, 'cached_full_gradient_points': len(records)})
    require(len(preflight) == 8, 'Eight arm/state preflights required')
    return cfg, range_module, pool, base, world, settings, starts, preflight, shared_cache, shared_records, lbfgs_direction, old_cfg


def projected_diagnostic(jac, gradient, active, modes):
    full = jac.reshape(64, -1).T; scales = np.linalg.norm(full, axis=0)
    require(np.all(scales > 1e-12), 'Zero current physical mode')
    normed = full/scales; matrix = normed[:, active]; eig, vectors = np.linalg.eigh(matrix.T@matrix)
    keep = eig > 1e-10*eig[-1]; require(int(keep.sum()) == len(active), 'Current arm metric rank deficient')
    inverse = (vectors[:, keep]/eig[keep])@vectors[:, keep].T
    flat = gradient.reshape(-1); projected = (matrix@inverse@matrix.T@flat).reshape(256, 3)
    residual = gradient-projected; loads = -normed.T@flat
    row = {'active_span_max_atom_kcal_mol_A': float(np.linalg.norm(projected, axis=1).max()),
        'active_span_RMS_cartesian_kcal_mol_A': float(np.sqrt(np.mean(projected**2))),
        'active_span_L2_kcal_mol_A': float(np.linalg.norm(projected)),
        'full_gradient_L2_kcal_mol_A': float(np.linalg.norm(gradient)),
        'outside_active_span_L2_kcal_mol_A': float(np.linalg.norm(residual)),
        'outside_active_span_max_atom_kcal_mol_A': float(np.linalg.norm(residual, axis=1).max()),
        'normalized_Gram_eigenvalues': eig.tolist(), 'rank': len(active), 'stationary_minimum_claimed': False,
        'mode_force_loads_kcal_mol_A': [{'mode_index': k, 'name': m['name'], 'family': m['family'],
            'active': k in active, 'force_load': float(loads[k])} for k, m in enumerate(modes)]}
    row['low_active_span_force'] = bool(row['active_span_max_atom_kcal_mol_A'] <= .1 and row['active_span_RMS_cartesian_kcal_mol_A'] <= .03)
    return row


def feasible_coordinates(start, q, positions, jac, gradient):
    """One downhill physical-coordinate candidate per mode; geometry only."""
    kin, guard, active = start['source']['kin'], start['guard'], start['active']
    gq = np.einsum('kij,ij->k', jac, gradient); candidates, records = [], []
    for k in sorted(active):
        speed = float(np.linalg.norm(jac[k], axis=1).max())
        record = {'mode_index': k, 'mode_name': kin.modes[k]['name'], 'coordinate_gradient': float(gq[k]),
            'tangent_max_atom_A_per_unit': speed, 'geometric_trials': [], 'admitted': False}
        records.append(record)
        if speed <= 1e-14 or gq[k] == 0:
            record['reason'] = 'zero_tangent_or_coordinate_gradient'; continue
        step = -float(np.sign(gq[k]))*.05/speed
        for multiplier in (1., .5, .25, .125):
            delta = multiplier*step; trial_q = q.copy(); trial_q[k] += delta
            trial_positions, _, trial_jac, _ = kin.evaluate(trial_q)
            check = guard.check(trial_q, trial_positions, positions)
            record['geometric_trials'].append({'multiplier': multiplier, 'coordinate_step': delta, 'geometry': check})
            if not check['pass']:
                continue
            slope = float(gq[k]*delta); require(slope < 0, 'Fallback coordinate is not own-energy descent')
            trial_v = np.linalg.solve(start['transform'], trial_q[active]-start['q'][active])
            require(np.allclose(start['q'][active]+start['transform']@trial_v, trial_q[active], rtol=0, atol=1e-12),
                    'Fallback physical/offset coordinates differ')
            record.update(admitted=True, multiplier=multiplier, coordinate_step=delta,
                predicted_own_decrease_kcal_mol=-slope, full_q=trial_q.tolist(), offset_v_A=trial_v.tolist())
            candidates.append({'record': record, 'q': trial_q, 'v': trial_v, 'positions': trial_positions,
                'jac': trial_jac, 'guard': check, 'slope': slope})
            break  # first admissible multiplier is the largest prescribed one
    candidates.sort(key=lambda c: (-c['record']['predicted_own_decrease_kcal_mol'], c['record']['mode_index']))
    return candidates, records


def search(base, start, arm, evaluator, out, lbfgs_direction):
    source, active, transform = start['source'], start['active'], start['transform']
    kin, guard = source['kin'], start['guard']; n = len(active)
    q0 = start['q'].copy(); q = q0.copy(); v = np.zeros(n)
    positions, jac, gradient = (start[key].copy() for key in ('x', 'jac', 'gradient'))
    point, point_pin = start['point'], start['point_pin']; origin = start['origin']; history = []
    gv = transform.T@np.einsum('kij,ij->k', jac[active], gradient)
    current_guard = start['start_guard']
    result = {'protocol_id': PROTOCOL, 'source_id': source['id'], 'metal': base.METAL, 'arm': arm,
        'status': 'unavailable', 'candidate_admitted': False, 'active_indices': active,
        'active_names': [kin.modes[k]['name'] for k in active], **start['metric'],
        'coordinate_offset_q64': q0.tolist(), 'coordinate_rule': 'q_active=q_start_active+T_arm@v; v initially zero',
        'inherited_seven_coordinate_scaling_T': start['search']['coordinate_scaling_T'],
        'authoritative_origin': start['search']['authoritative_origin'], 'original_point': start['origin_pin'],
        'previous_search': start['search_pin'], 'continuation_start': start['point_pin'],
        'origin_guard': start['origin_guard'], 'continuation_start_guard': current_guard,
        'initial_history_length': 0, 'initial_new_evaluations': 0, 'initial_offset_v_A': v.tolist(),
        'source_rebased': False, 'trace': [], 'accepted_steps': 0, 'accepted_fallback_steps': 0,
        'fallback_rule': FALLBACK, 'stationary_minimum_claimed': False}
    try:
        while True:
            diagnostic = projected_diagnostic(jac, gradient, active, kin.modes)
            row = {'accepted_step': result['accepted_steps'], 'offset_v_A': v.tolist(), 'full_q': q.tolist(), 'point': point_pin,
                'work_from_origin': base.relative(point['components'], origin['components']),
                'work_from_continuation_start': base.relative(point['components'], start['point']['components']),
                'geometry': current_guard, 'force_diagnostic': diagnostic}
            result['trace'].append(row)
            with (out/'TRACE.jsonl').open('a') as handle:
                handle.write(json.dumps(row, sort_keys=True, allow_nan=False)+'\n')
            if diagnostic['low_active_span_force']:
                termination = 'low_active_span_force'; break
            if result['accepted_steps'] >= 12:
                termination = 'accepted_step_limit'; break
            if evaluator.cost['new_composite_evaluations'] >= 24:
                termination = 'evaluation_limit'; break
            direction = lbfgs_direction(gv, history)
            if float(gv@direction) >= -1e-14:
                history = []; direction = -gv
            tangent = np.einsum('k,kij->ij', transform@direction, jac[active])
            speed = float(np.linalg.norm(tangent, axis=1).max())
            if speed <= 1e-14:
                termination = 'zero_active_direction'; break
            direction *= min(1., .05/speed); slope = float(gv@direction)
            require(slope < 0, 'Own-energy direction is not descent')
            row['trials'] = []; accepted = None; used_fallback = False
            for multiplier in base.SETTINGS['trial_multipliers']:
                trial_v = v+multiplier*direction; trial_q = q0.copy(); trial_q[active] = q0[active]+transform@trial_v
                trial_positions, _, trial_jac, _ = kin.evaluate(trial_q)
                check = guard.check(trial_q, trial_positions, positions)
                trial = {'multiplier': multiplier, 'offset_v_A': trial_v.tolist(), 'full_q': trial_q.tolist(),
                    'geometry': check, 'molecular_evaluation_requested': False, 'accepted': False}
                row['trials'].append(trial)
                if not check['pass']:
                    trial['reason'] = 'geometry_or_bound'; continue
                if evaluator.cost['new_composite_evaluations'] >= 24 and evaluator.key(trial_positions) not in evaluator.cache:
                    trial['reason'] = 'evaluation_limit'; break
                trial['molecular_evaluation_requested'] = True  # cache-energy trials also prohibit fallback entry
                next_point, next_pin = evaluator.evaluate(trial_positions)
                change = base.relative(next_point['components'], point['components'])
                trial.update(point=next_pin, work_from_current=change, armijo_rhs_kcal_mol=1e-4*multiplier*slope)
                if change['composite_kcal_mol'] <= trial['armijo_rhs_kcal_mol']:
                    trial['accepted'] = True
                    accepted = (trial_v, trial_q, trial_positions, trial_jac, next_point, next_pin, check); break
                trial['reason'] = 'own_composite_Armijo'
            fallback_allowed = (accepted is None and len(row['trials']) == 4
                and all(t.get('reason') == 'geometry_or_bound' and not t['molecular_evaluation_requested'] for t in row['trials']))
            row['fallback'] = {'entered': fallback_allowed, 'evaluated_candidates': []}
            if fallback_allowed:
                candidates, proposals = feasible_coordinates(start, q, positions, jac, gradient)
                row['fallback'].update(coordinate_proposals=proposals,
                    ranked_mode_indices=[c['record']['mode_index'] for c in candidates])
                for candidate in candidates[:2]:
                    if evaluator.cost['new_composite_evaluations'] >= 24 and evaluator.key(candidate['positions']) not in evaluator.cache:
                        row['fallback']['stop'] = 'evaluation_limit'; break
                    trial = {'mode_index': candidate['record']['mode_index'], 'accepted': False,
                        'predicted_own_decrease_kcal_mol': -candidate['slope'],
                        'armijo_rhs_kcal_mol': 1e-4*candidate['slope'], 'molecular_evaluation_requested': True}
                    row['fallback']['evaluated_candidates'].append(trial)
                    next_point, next_pin = evaluator.evaluate(candidate['positions'])
                    change = base.relative(next_point['components'], point['components'])
                    trial.update(point=next_pin, work_from_current=change)
                    if change['composite_kcal_mol'] <= trial['armijo_rhs_kcal_mol']:
                        trial['accepted'] = True; used_fallback = True
                        accepted = (candidate['v'], candidate['q'], candidate['positions'], candidate['jac'],
                                    next_point, next_pin, candidate['guard']); break
                    trial['reason'] = 'own_composite_Armijo'
            write(out/('ITERATION_%02d.json' % result['accepted_steps']), row)
            if accepted is None:
                termination = 'evaluation_limit' if evaluator.cost['new_composite_evaluations'] >= 24 else 'no_admissible_decreasing_trial'; break
            new_v, q, positions, jac, point, point_pin, current_guard = accepted
            gradient = np.load(verify(point['gradient_kcal_mol_A']), allow_pickle=False)
            new_gv = transform.T@np.einsum('kij,ij->k', jac[active], gradient)
            if used_fallback:
                history = []; result['accepted_fallback_steps'] += 1
            else:
                displacement, difference = new_v-v, new_gv-gv; sy = float(displacement@difference)
                if sy > 1e-10*np.linalg.norm(displacement)*np.linalg.norm(difference) and sy > 1e-16:
                    history.append((displacement, difference, 1./sy)); history = history[-8:]
            v, gv = new_v, new_gv; result['accepted_steps'] += 1
        final_check = guard.check(q, positions, positions); require(final_check['pass'], 'Terminal physical gate failed')
        work = base.relative(point['components'], origin['components'])
        incremental = base.relative(point['components'], start['point']['components'])
        require(work['composite_kcal_mol'] <= 1e-8 and incremental['composite_kcal_mol'] <= 1e-8, 'Own terminal energy increased')
        result.update(status='finite_proposal_available', candidate_admitted=True, termination=termination, candidate=point_pin,
            full_q=q.tolist(), offset_v_A=v.tolist(), final_geometry=final_check, final_force_diagnostic=diagnostic,
            work_from_origin=work, work_from_continuation_start=incremental, preserved_boundary_flags=final_check['boundary_flags'])
    except Exception as exc:
        result.update(status='unavailable', candidate_admitted=False, termination='backend_or_numerical_failure',
            failure=str(exc), traceback=traceback.format_exc(), last_archived_valid_point=point_pin)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('preflight', 'scout', 'remaining'), required=True)
    for name in ('config', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args(); require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic()
    result = {'protocol_id': PROTOCOL, 'stage': args.stage, 'status': 'unavailable',
        'implementation': pin(__file__), 'config': pin(args.config), 'preflight': {'status': 'not_passed'},
        'rows': [], 'points': [], 'new_initialization_molecular_calls': 0, 'source_rebased': False,
        'pool_comparison_performed': False, 'stationary_minimum_claimed': False, 'fallback_rule': FALLBACK,
        'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname()}
    write(out/'START.json', result)
    try:
        (cfg, range_module, pool, base, world, settings, starts, checks, shared_cache,
         shared_records, lbfgs_direction, old_cfg) = inputs(args.config)
        range_scout = read(cfg['range_scout_result'])
        result.update(plan=cfg['plan'], settings=settings, range_scout_result=cfg['range_scout_result'],
            range_remaining_result=cfg['range_remaining_result'], old_pool_result=range_scout['old_pool_result'],
            frozen_scout_implementation=range_scout['frozen_scout_implementation'], frozen_range_implementation=range_scout['implementation'],
            profile_result=old_cfg['profile_result'], fd_result=old_cfg['fd_result'], helpers=world['profile']['helpers'],
            native_helper=range_scout['native_helper'], solvent_helper=range_scout['solvent_helper'],
            underlying_molecular_protocol=base.PROTOCOL,
            preflight={'status': 'passed', 'declared_arm_states': 8, 'rows': checks, 'molecular_calls': 0, 'model_loaded': False})
        write(out/'PREFLIGHT.json', result['preflight']); result['preflight_result'] = pin(out/'PREFLIGHT.json')
        if args.stage == 'preflight':
            result['status'] = 'complete'
        else:
            require(int(os.environ['SLURM_CPUS_PER_TASK']) >= 2
                    and os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'), 'Allocated GPU/two solvent workers required')
            if args.stage == 'remaining':
                scout = read(cfg['basis_scout_result']); scout_cfg = read(scout['config'])
                require(scout['protocol_id'] == PROTOCOL and scout['stage'] == 'scout' and scout['status'] == 'complete'
                        and scout['continuation_execution_and_physical_admission'] and scout['plan'] == cfg['plan']
                        and scout['settings'] == settings and scout['fallback_rule'] == FALLBACK
                        and scout['implementation'] == scout_cfg['implementation']
                        and scout['implementation']['sha256'] == cfg['implementation']['sha256']
                        and all(scout[key] == cfg[key] for key in ('range_scout_result', 'range_remaining_result'))
                        and len(scout['rows']) == 2
                        and [(r['source_id'], r['metal'], r['arm']) for r in scout['rows']] == [('7CCO', 'La', a) for a in ARMS],
                        'Both matched scout arms required')
                for row in scout['rows']:
                    require(row['candidate_admitted'] and row['status'] == 'complete'
                            and read(row['result']) == {k: v for k, v in row.items() if k != 'result'}, 'Basis scout row gate failed')
                    proposal = read(row['search']); require(proposal['candidate_admitted'] and proposal['final_geometry']['pass'], 'Basis scout physical gate failed')
                result['basis_scout_result'] = cfg['basis_scout_result']
            states = STATES[:1] if args.stage == 'scout' else STATES[1:]
            tasks = [(s, m, arm) for s, m in states for arm in ARMS]
            result['declared_tasks'] = [{'source_id': s, 'metal': m, 'arm': a} for s, m, a in tasks]
            evaluator_type = range_module.evaluator_class(base); warm = None
            for source_id, metal, arm in tasks:
                directory = out/(source_id+'_'+metal+'_'+arm); directory.mkdir()
                base.SOURCE, base.METAL = source_id, metal; state = (source_id, metal)
                start = dict(starts[(source_id, metal, arm)], cache=shared_cache[state], records=shared_records[state])
                evaluator = None; prior_calls = warm.calls if warm is not None else 0
                row = {'source_id': source_id, 'metal': metal, 'arm': arm, 'status': 'unavailable', 'candidate_admitted': False,
                    'physical_multiplicity': 1 if metal == 'La' else 6, 'GFN2_effective_multiplicity': 1,
                    'previous_search': start['search_pin'], 'continuation_start': start['point_pin']}
                try:
                    evaluator = evaluator_type(world, start['source'], directory, start); evaluator.native = warm
                    proposal = search(base, start, arm, evaluator, directory, lbfgs_direction)
                    write(directory/'SEARCH.json', proposal)
                    row.update(status='complete' if proposal['candidate_admitted'] else 'unavailable',
                        candidate_admitted=proposal['candidate_admitted'], search=pin(directory/'SEARCH.json'))
                except Exception as exc:
                    row.update(failure=str(exc), traceback=traceback.format_exc())
                if evaluator is not None:
                    warm = evaluator.native
                    # Search starts stay fixed; only exact compatible molecular evaluations are shared.
                    shared_cache[state] = dict(evaluator.cache); shared_records[state] = list(evaluator.records)
                    row.update(actual_cost=evaluator.cost, points=evaluator.records, requests=evaluator.requests)
                    result['points'].extend(evaluator.records)
                    if warm is not None:
                        row['native_runtime'] = {'versions': warm.versions, 'device': warm.device,
                            'worker_prior_forwards': prior_calls, 'worker_cumulative_forwards': warm.calls,
                            'actual_model_forwards_this_search': warm.calls-prior_calls, 'torch_threads': warm.torch.get_num_threads()}
                        require(warm.calls-prior_calls == evaluator.cost['MACE_calls_attempted'], 'Per-arm native call count differs')
                write(directory/'RESULT.json', row); row['result'] = pin(directory/'RESULT.json'); result['rows'].append(row)
                print(json.dumps({'source_id': source_id, 'metal': metal, 'arm': arm, 'status': row['status']}), flush=True)
            result['status'] = 'complete' if all(r['candidate_admitted'] for r in result['rows']) else 'unavailable'
            result['continuation_execution_and_physical_admission'] = result['status'] == 'complete'
            if warm is not None:
                result['native_runtime'] = {'actual_model_forwards': warm.calls, 'model_loaded_once': True,
                    'peak_cuda_allocated_bytes': warm.torch.cuda.max_memory_allocated()}
        verify(result['config'])
    except Exception as exc:
        result.update(status='execution_failure', failure=str(exc), traceback=traceback.format_exc())
    result['actual_cost'] = {key: sum(r.get('actual_cost', {}).get(key, 0) for r in result['rows']) for key in
        ('new_composite_evaluations', 'MACE_calls_attempted', 'GFN2_calls_attempted', 'MACE_call_seconds',
         'GFN2_call_seconds', 'model_load_seconds', 'cache_hits')}
    limit = 2 if args.stage == 'scout' else 6 if args.stage == 'remaining' else 0
    require(result['actual_cost']['new_composite_evaluations'] <= 24*limit and result['actual_cost']['MACE_calls_attempted'] <= 24*limit
            and result['actual_cost']['GFN2_calls_attempted'] <= 48*limit, 'Matched arm call ceiling exceeded')
    result['elapsed_seconds'] = time.monotonic()-tick
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds']*int(os.environ['SLURM_CPUS_PER_TASK'])
    result['peak_host_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['allocation'] = {k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_CPU', 'CUDA_VISIBLE_DEVICES')}
    write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}), flush=True)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
