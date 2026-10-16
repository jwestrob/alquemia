#!/usr/bin/env python3
"""One declared range-only continuation, using exact old terminal gradients.

CLI: --stage preflight|scout|remaining --config CONFIG --output FRESH
Config pins: implementation, plan, old_pool_result. Remaining additionally pins
range_scout_result. Every stage preflights all four actual old starts without
loading a model. Numerical execution is Slurm-only. No pool scoring is done here.
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

REPO = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
PROTOCOL = 'lbt3_whole_composite_range04_continuation_v1'
SCOUT_SHA = 'ee1db81bdf3db90c303f824d7a55b52b00c2bb64f2422903aed60c25931d8ca1'
POOL_SHA = '505ae7e59ba6abcd2575536f40291bb1c7932a911e21e3aac6268d865885ae54'
TASKS = (('7CCO', 'La'), ('7CCO', 'Dy'), ('7CCN', 'La'), ('7CCN', 'Dy'))


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


def range_guard(base, source, active, Geometry, geometry_settings, metal):
    """Reuse every old physical check; replace only its hardcoded source extent.

Public factory for terminal pool admission. Caller initializes base.np and this
module's np. The original source is never replaced by a continuation endpoint.
"""
    require(base.SETTINGS['all_atom_displacement_A'] == .20 and metal in ('La', 'Dy'), 'Unexpected parent guard')

    class RangeGuard(base.GeometryGuard):
        def check(self, q, positions, previous):
            row = super().check(q, positions, previous)
            row['old_020_guard_pass'] = row['pass']
            eps = base.SETTINGS['displacement_roundoff_A']
            row['pass'] = bool(row['peptide_guard']['pass'] and not row['all256_new_worsened_severe_contacts']
                and row['maximum_source_displacement_A'] <= .40 + eps
                and row['target_displacement_A'] <= .20 + eps
                and row['maximum_trial_displacement_A'] <= .05 + eps
                and row['maximum_torsion_radian'] <= .20 + 1e-12
                and row['fixed_atom_max_change_A'] <= 1e-12
                and row['peptide_guard']['source_bond_max_error_A'] <= 1e-9
                and row['bond_angle_cosine_max_change'] <= 1e-10
                and row['protected_group_distance_max_change_A'] <= 1e-9)
            row['boundary_flags']['source_displacement'] = row['maximum_source_displacement_A'] >= .40 - 1e-5
            row['source_displacement_bound_A'] = .40
            row['source_baseline'] = source['prepared']['physical_atoms']
            return row

    previous_metal = base.METAL
    try:
        base.METAL = metal
        return RangeGuard(source, active, Geometry, geometry_settings)
    finally:
        base.METAL = previous_metal


def saved_point(pool, world, source, metal, record):
    point = read(record)
    x, _ = pool.point_components(world, source, metal, point)
    gradient = np.load(verify(point['gradient_kcal_mol_A']), allow_pickle=False)
    require(gradient.shape == (256, 3) and np.isfinite(gradient).all(), 'Saved gradient invalid')
    if 'reused_origin' in point:
        expected = np.asarray(world['corrections'][(source['id'], metal)]['composite_gradient']['kcal_mol_A'])
    else:
        mace = read(point['MACE'])
        vacuum, water = [read(point['GFN2'][medium]) for medium in ('vacuum', 'ALPB(Water)')]
        expected = (-np.asarray(mace['forces_eV_A']) * pool.EV_TO_KCAL
                    + np.asarray(water['gradient']['kcal_mol_A']) - np.asarray(vacuum['gradient']['kcal_mol_A']))
    require(np.allclose(gradient, expected, rtol=0, atol=1e-10), 'Saved component/full-gradient replay differs')
    return point, x, gradient


def inputs(config_path):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    cfg = json.loads(config_path.read_text()); verify(cfg['plan']); verify(cfg['implementation'])
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Range adapter differs')
    old = read(cfg['old_pool_result']); old_cfg = read(old['config'])
    require(old['status'] == 'complete' and old['stage'] == 'run'
            and old['implementation'] == old_cfg['implementation']
            and old['implementation']['sha256'] == POOL_SHA
            and old['design'] == dict({k: old_cfg[k] for k in ('plan', 'profile_result', 'fd_result', 'scout_result', 'remaining_result')},
                                     implementation_sha256=POOL_SHA), 'Old complete pool/config lineage differs')
    pool = module('lbt3_range_old_pool', old['implementation'])
    inherited_cfg, world = pool.inputs(verify(old['config']))
    candidates, aliases, cells, slots, missing = pool.assemble(world)
    require(all(s['status'] == 'admitted' for s in slots) and old['searches'] == slots
            and old['aliases'] == aliases and old['candidates'] == [{k: v for k, v in c.items() if k != 'x'} for c in candidates]
            and old['required_cross_cells'] == [{'candidate_id': i, 'metal': m} for i, m in missing], 'Old pool candidate/row replay differs')
    for row in old['cross_rows']:
        require(row['status'] == 'accepted' and read(row['cell']) == {k: v for k, v in row.items() if k != 'cell'},
                'Old cross-cell receipt differs')
    require(len(old['cross_rows']) == len(missing)
            and {(r['candidate_id'], r['metal']) for r in old['cross_rows']} == set(missing), 'Old paired pool incomplete')
    scout = world['parents']['scout']
    require(scout['implementation']['sha256'] == SCOUT_SHA, 'Frozen search implementation differs')
    base = module('lbt3_range_frozen_search', scout['implementation'])
    global np
    import numpy as np
    base.np = np
    require(base.SETTINGS == scout['settings'], 'Original search settings differ')
    settings = dict(base.SETTINGS, all_atom_displacement_A=.40)
    active = world['profile']['active_candidate_basis']['mode_indices_selection_order']
    require(world['profile']['active_candidate_basis']['mode_names_selection_order'] == base.NAMES, 'Active basis differs')
    solvent = read(world['profile']['solvent_result'])
    world.update(modes=read(world['profile']['modes_result']),
        solvent_module=module('lbt3_range_solvent', solvent['implementation']), fd_module=world['fd'],
        corrections={(r['source_id'], r['metal']): r for r in solvent['corrections']})
    require(pin(sys.executable)['sha256'] == world['mcfg']['python']['sha256'], 'Native Python differs')
    from scaffold_accommodation import Geometry, SETTINGS as GEOMETRY_SETTINGS, lbfgs_direction
    starts = {}; checks = []
    for slot in sorted(world['slots'], key=lambda r: TASKS.index((r['source_id'], r['metal']))):
        source_id, metal = slot['source_id'], slot['metal']; source = world['sources'][source_id]
        search = read(slot['search']); cache = {}; records = []
        for record in slot['points']:
            point, x, gradient = saved_point(pool, world, source, metal, record)
            if point['cache_key'] in cache:
                require(read(cache[point['cache_key']])['components'] == point['components'], 'Conflicting old point cache')
            else:
                cache[point['cache_key']] = record; records.append(record)
        point, x, gradient = saved_point(pool, world, source, metal, search['candidate'])
        origin_pin = search['trace'][0]['point']; origin, origin_x, _ = saved_point(pool, world, source, metal, origin_pin)
        require(np.array_equal(origin_x, source['x']) and search['candidate'] in records and origin_pin in records,
                'Exact original/start point absent from cache')
        q, u, transform = (np.asarray(search[key], dtype=float) for key in ('full_q', 'u_A', 'coordinate_scaling_T'))
        require(q.shape == (64,) and u.shape == (7,) and transform.shape == (7, 7)
                and all(np.isfinite(a).all() for a in (q, u, transform))
                and np.count_nonzero(q[[k for k in range(64) if k not in active]]) == 0
                and np.allclose(q[active], transform @ u, rtol=0, atol=1e-12), 'Original coordinates/whitening changed')
        replay, _, jac, _ = source['kin'].evaluate(q)
        zero, _, jac0, _ = source['kin'].evaluate(np.zeros(64))
        matrix = jac0[active].reshape(7, -1).T
        require(np.allclose(replay, x, rtol=0, atol=1e-12) and np.allclose(zero, source['x'], rtol=0, atol=1e-12)
                and np.allclose((matrix @ transform).T @ (matrix @ transform), np.eye(7), rtol=0, atol=1e-10)
                and np.allclose(np.linalg.norm(matrix, axis=0), search['source_column_norms'], rtol=0, atol=1e-12),
                'Original metric/terminal kinematics replay differs')
        guard = range_guard(base, source, active, Geometry, GEOMETRY_SETTINGS, metal)
        origin_guard = guard.check(np.zeros(64), source['x'], source['x']); start_guard = guard.check(q, x, x)
        require(origin_guard['pass'] and start_guard['pass'] and start_guard['old_020_guard_pass'], 'Old origin/start physical gate failed')
        checks.append({'source_id': source_id, 'metal': metal, 'pass': True, 'previous_search': slot['search'],
            'continuation_start': search['candidate'], 'original_point': origin_pin, 'cached_full_gradient_points': len(records),
            'original_whitening_reused': True, 'initial_history_length': 0, 'initial_new_evaluations': 0,
            'gradient_replay_max_error_A_units': 'component gradient replay absolute tolerance 1e-10 kcal/mol/angstrom',
            'start_kinematic_replay_max_error_A': float(np.max(np.abs(replay-x))),
            'original_guard': origin_guard, 'start_guard': start_guard})
        starts[(source_id, metal)] = dict(source=source, search=search, search_pin=slot['search'], point=point,
            point_pin=search['candidate'], origin=origin, origin_pin=origin_pin, x=x, gradient=gradient,
            q=q, u=u, transform=transform, jac=jac, cache=cache, records=records,
            guard=guard, origin_guard=origin_guard, start_guard=start_guard)
    return cfg, old, inherited_cfg, pool, base, world, active, settings, starts, checks, lbfgs_direction


def evaluator_class(base):
    class ContinuationEvaluator(base.CompositeEvaluator):
        def __init__(self, world, source, out, start):
            self.start = start
            super().__init__(world, source, out)

        def save_origin(self):
            # Parent initializer calls here before creating an origin key/file.
            self.recipe['physical_multiplicity'] = self.base['multiplicity']
            require(self.base['multiplicity'] == self.base['spin_multiplicity'] == (1 if base.METAL == 'La' else 6)
                    and self.start['origin']['recipe'] == self.recipe
                    and self.start['point']['recipe'] == self.recipe, 'Scoped physical cache/state differs')
            self.cache.update(self.start['cache']); self.records.extend(self.start['records'])
            require(self.key(self.start['x']) == self.start['point']['cache_key'], 'Warm start cache key differs')
            return self.start['origin'], self.start['origin_pin']

    return ContinuationEvaluator


def search(base, start, active, evaluator, out, lbfgs_direction):
    """Frozen own-energy loop, with saved q/u/metric/gradient initialization."""
    source = start['source']; kin = source['kin']; transform = start['transform']
    q, u, positions, jac = [start[key].copy() for key in ('q', 'u', 'x', 'jac')]
    gradient = start['gradient'].copy(); guard = start['guard']; history = []
    point, point_pin = start['point'], start['point_pin']; origin = start['origin']
    gu = transform.T @ np.einsum('kij,ij->k', jac[active], gradient)
    current_guard = start['start_guard']; old = start['search']
    result = {'protocol_id': PROTOCOL, 'source_id': source['id'], 'metal': base.METAL,
        'status': 'unavailable', 'candidate_admitted': False, 'active_indices': active,
        'active_names': base.NAMES, 'coordinate_scaling_T': old['coordinate_scaling_T'],
        'source_column_norms': old['source_column_norms'],
        'source_normalized_Gram_eigenvalues': old['source_normalized_Gram_eigenvalues'],
        'q0_kinematic_replay_max_error_A': old['q0_kinematic_replay_max_error_A'],
        'authoritative_origin': old['authoritative_origin'], 'original_point': start['origin_pin'],
        'previous_search': start['search_pin'], 'continuation_start': start['point_pin'],
        'origin_guard': start['origin_guard'], 'continuation_start_guard': current_guard,
        'initial_history_length': 0, 'initial_new_evaluations': 0, 'trace': [], 'accepted_steps': 0,
        'stationary_minimum_claimed': False, 'source_rebased': False,
        'origin_coordinate_policy': 'exact original source and old terminal point; saved original whitening, no new initialization calls'}
    try:
        while True:
            diagnostic = base.projected_diagnostic(jac, gradient, active, kin.modes)
            row = {'accepted_step': result['accepted_steps'], 'u_A': u.tolist(), 'full_q': q.tolist(),
                'point': point_pin, 'work_from_origin': base.relative(point['components'], origin['components']),
                'work_from_continuation_start': base.relative(point['components'], start['point']['components']),
                'geometry': current_guard, 'force_diagnostic': diagnostic}
            result['trace'].append(row)
            with (out / 'TRACE.jsonl').open('a') as handle:
                handle.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')
            if diagnostic['low_active_span_force']:
                termination = 'low_active_span_force'; break
            if result['accepted_steps'] >= 12:
                termination = 'accepted_step_limit'; break
            if evaluator.cost['new_composite_evaluations'] >= 24:
                termination = 'evaluation_limit'; break
            direction = lbfgs_direction(gu, history)
            if float(gu @ direction) >= -1e-14:
                history = []; direction = -gu
            physical_tangent = np.einsum('k,kij->ij', transform @ direction, jac[active])
            speed = float(np.linalg.norm(physical_tangent, axis=1).max())
            if speed <= 1e-14:
                termination = 'zero_active_direction'; break
            direction *= min(1., .05 / speed)
            slope = float(gu @ direction); require(slope < 0, 'Own-energy direction is not descent')
            accepted = None; row['trials'] = []
            for multiplier in base.SETTINGS['trial_multipliers']:
                trial_u = u + multiplier * direction; trial_q = np.zeros(64); trial_q[active] = transform @ trial_u
                trial_positions, _, trial_jac, _ = kin.evaluate(trial_q)
                check = guard.check(trial_q, trial_positions, positions)
                trial = {'multiplier': multiplier, 'u_A': trial_u.tolist(), 'full_q': trial_q.tolist(),
                    'geometry': check, 'molecular_evaluation_requested': False, 'accepted': False}
                row['trials'].append(trial)
                if not check['pass']:
                    trial['reason'] = 'geometry_or_bound'; continue
                if evaluator.cost['new_composite_evaluations'] >= 24 and evaluator.key(trial_positions) not in evaluator.cache:
                    trial['reason'] = 'evaluation_limit'; break
                trial['molecular_evaluation_requested'] = True
                next_point, next_pin = evaluator.evaluate(trial_positions)
                change = base.relative(next_point['components'], point['components'])
                trial.update(point=next_pin, work_from_current=change, armijo_rhs_kcal_mol=1e-4 * multiplier * slope)
                if change['composite_kcal_mol'] <= trial['armijo_rhs_kcal_mol']:
                    trial['accepted'] = True
                    accepted = (trial_u, trial_q, trial_positions, trial_jac, next_point, next_pin, check)
                    break
                trial['reason'] = 'own_composite_Armijo'
            write(out / ('ITERATION_%02d.json' % result['accepted_steps']), row)
            if accepted is None:
                termination = 'evaluation_limit' if evaluator.cost['new_composite_evaluations'] >= 24 else 'no_admissible_decreasing_trial'
                break
            new_u, q, positions, jac, point, point_pin, current_guard = accepted
            gradient = np.load(verify(point['gradient_kcal_mol_A']), allow_pickle=False)
            new_gu = transform.T @ np.einsum('kij,ij->k', jac[active], gradient)
            displacement, difference = new_u - u, new_gu - gu
            sy = float(displacement @ difference)
            if sy > 1e-10 * np.linalg.norm(displacement) * np.linalg.norm(difference) and sy > 1e-16:
                history.append((displacement, difference, 1. / sy)); history = history[-8:]
            u, gu = new_u, new_gu; result['accepted_steps'] += 1
        final_check = guard.check(q, positions, positions)
        require(final_check['pass'], 'Last accepted geometry failed final guard')
        work = base.relative(point['components'], origin['components'])
        incremental = base.relative(point['components'], start['point']['components'])
        require(work['composite_kcal_mol'] <= 1e-8 and incremental['composite_kcal_mol'] <= 1e-8,
                'Last accepted own energy increased')
        result.update(status='finite_proposal_available', candidate_admitted=True, termination=termination,
            candidate=point_pin, full_q=q.tolist(), u_A=u.tolist(), final_geometry=final_check,
            work_from_origin=work, work_from_continuation_start=incremental,
            final_force_diagnostic=diagnostic, preserved_boundary_flags=final_check['boundary_flags'])
    except Exception as exc:
        result.update(status='unavailable', candidate_admitted=False, termination='backend_or_numerical_failure',
            failure=str(exc), traceback=traceback.format_exc(), last_archived_valid_point=point_pin)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('preflight', 'scout', 'remaining'), required=True)
    for name in ('config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic()
    result = {'protocol_id': PROTOCOL, 'stage': args.stage, 'status': 'unavailable',
        'implementation': pin(__file__), 'config': pin(args.config), 'preflight': {'status': 'not_passed'},
        'rows': [], 'points': [], 'new_initialization_molecular_calls': 0, 'source_rebased': False,
        'pool_comparison_performed': False, 'stationary_minimum_claimed': False,
        'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname()}
    write(out / 'START.json', result)
    try:
        cfg, old, inherited_cfg, pool, base, world, active, settings, starts, checks, lbfgs_direction = inputs(args.config)
        result.update(plan=cfg['plan'], old_pool_result=cfg['old_pool_result'], settings=settings,
            underlying_molecular_protocol=base.PROTOCOL, frozen_scout_implementation=world['parents']['scout']['implementation'],
            frozen_pool_implementation=old['implementation'], profile_result=inherited_cfg['profile_result'],
            fd_result=inherited_cfg['fd_result'], native_helper=world['parents']['scout']['native_helper'],
            solvent_helper=world['parents']['scout']['solvent_helper'], helpers=world['profile']['helpers'],
            preflight={'status': 'passed', 'declared_states': 4, 'rows': checks, 'molecular_calls': 0, 'model_loaded': False})
        write(out / 'PREFLIGHT.json', result['preflight']); result['preflight_result'] = pin(out / 'PREFLIGHT.json')
        if args.stage == 'preflight':
            result['status'] = 'complete'
        else:
            require(int(os.environ['SLURM_CPUS_PER_TASK']) >= 2
                    and os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'), 'Allocated GPU/two solvent workers required')
            if args.stage == 'remaining':
                scout = read(cfg['range_scout_result']); scout_cfg = read(scout['config'])
                require(scout['protocol_id'] == PROTOCOL and scout['stage'] == 'scout' and scout['status'] == 'complete'
                        and scout['continuation_execution_and_physical_admission']
                        and scout['plan'] == cfg['plan'] and scout['old_pool_result'] == cfg['old_pool_result']
                        and scout['settings'] == settings and scout['implementation'] == scout_cfg['implementation']
                        and scout['implementation']['sha256'] == cfg['implementation']['sha256']
                        and len(scout['rows']) == 1 and scout['rows'][0]['source_id'] == '7CCO'
                        and scout['rows'][0]['metal'] == 'La' and scout['rows'][0]['candidate_admitted'],
                        'Range scout execution/physical gate failed')
                admitted = read(scout['rows'][0]['search'])
                require(admitted['candidate_admitted'] and admitted['final_geometry']['pass'], 'Range scout terminal unavailable')
                result['range_scout_result'] = cfg['range_scout_result']
            tasks = TASKS[:1] if args.stage == 'scout' else TASKS[1:]
            result['declared_tasks'] = [{'source_id': s, 'metal': m} for s, m in tasks]
            evaluator_type = evaluator_class(base); warm = None
            for source_id, metal in tasks:
                directory = out / (source_id + '_' + metal); directory.mkdir()
                base.SOURCE, base.METAL = source_id, metal
                start = starts[(source_id, metal)]; evaluator = None
                prior_calls = warm.calls if warm is not None else 0
                row = {'source_id': source_id, 'metal': metal, 'status': 'unavailable', 'candidate_admitted': False,
                    'physical_multiplicity': 1 if metal == 'La' else 6, 'GFN2_effective_multiplicity': 1,
                    'previous_search': start['search_pin'], 'continuation_start': start['point_pin']}
                try:
                    evaluator = evaluator_type(world, start['source'], directory, start); evaluator.native = warm
                    proposal = search(base, start, active, evaluator, directory, lbfgs_direction)
                    write(directory / 'SEARCH.json', proposal)
                    row.update(status='complete' if proposal['candidate_admitted'] else 'unavailable',
                        candidate_admitted=proposal['candidate_admitted'], search=pin(directory / 'SEARCH.json'))
                except Exception as exc:
                    row.update(failure=str(exc), traceback=traceback.format_exc())
                if evaluator is not None:
                    warm = evaluator.native
                    row.update(actual_cost=evaluator.cost, points=evaluator.records, requests=evaluator.requests)
                    result['points'].extend(evaluator.records)
                    if warm is not None:
                        row['native_runtime'] = {'versions': warm.versions, 'device': warm.device,
                            'worker_prior_forwards': prior_calls, 'worker_cumulative_forwards': warm.calls,
                            'actual_model_forwards_this_search': warm.calls-prior_calls,
                            'torch_threads': warm.torch.get_num_threads()}
                        require(warm.calls-prior_calls == evaluator.cost['MACE_calls_attempted'], 'Per-search native count differs')
                write(directory / 'RESULT.json', row); row['result'] = pin(directory / 'RESULT.json'); result['rows'].append(row)
                print(json.dumps({'source_id': source_id, 'metal': metal, 'status': row['status']}), flush=True)
            result['status'] = 'complete' if all(r['candidate_admitted'] for r in result['rows']) else 'unavailable'
            result['continuation_execution_and_physical_admission'] = result['status'] == 'complete'
            result['actual_cost'] = {key: sum(r.get('actual_cost', {}).get(key, 0) for r in result['rows']) for key in
                ('new_composite_evaluations', 'MACE_calls_attempted', 'GFN2_calls_attempted', 'MACE_call_seconds',
                 'GFN2_call_seconds', 'model_load_seconds', 'cache_hits')}
            require(result['actual_cost']['new_composite_evaluations'] <= 24*len(tasks)
                    and result['actual_cost']['MACE_calls_attempted'] <= 24*len(tasks)
                    and result['actual_cost']['GFN2_calls_attempted'] <= 48*len(tasks), 'Continuation call ceiling exceeded')
            if warm is not None:
                result['native_runtime'] = {'actual_model_forwards': warm.calls, 'model_loaded_once': True,
                    'peak_cuda_allocated_bytes': warm.torch.cuda.max_memory_allocated()}
        verify(result['config'])
    except Exception as exc:
        result.update(status='execution_failure', failure=str(exc), traceback=traceback.format_exc())
    result.setdefault('actual_cost', {'new_composite_evaluations': 0, 'MACE_calls_attempted': 0, 'GFN2_calls_attempted': 0})
    result['elapsed_seconds'] = time.monotonic()-tick
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds']*int(os.environ['SLURM_CPUS_PER_TASK'])
    result['peak_host_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['allocation'] = {k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_CPU', 'CUDA_VISIBLE_DEVICES')}
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}), flush=True)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
