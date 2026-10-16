#!/usr/bin/env python3
"""One bounded whole-peptide own-metal composite response scout; Slurm only.

CLI: python LBT3_ACCOMMODATE.py --stage scout --config CONFIG --output FRESH
Optional --repo-root overrides the explicit repository path below.
Config pins: plan, implementation, profile_result, fd_result. This first version
runs only 7CCO/La. Model, solvent runtime and helpers are inherited from the
accepted profile12512, derivative12521 and their exact origin receipts.
No repeated origin calls, added force-field energy, rescue or pool inference.
"""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import socket
import sys
import threading
import time
import traceback

PROTOCOL = 'lbt3_whole_own_composite_bounded_accommodation_v1'
REPO = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
SOURCE, METAL = '7CCO', 'La'
EV_TO_KCAL, HA_TO_KCAL = 23.06054783061903, 627.509474
SETTINGS = {'torsion_bound_radian': 0.20, 'target_displacement_A': 0.20,
    'all_atom_displacement_A': 0.20, 'all_atom_trial_step_A': 0.05,
    'maximum_accepted_steps': 12, 'maximum_new_composite_evaluations': 24,
    'history': 8, 'armijo': 1e-4, 'trial_multipliers': [1., .5, .25, .125],
    'metric_relative_cutoff': 1e-10, 'active_max_atom_kcal_mol_A': .1,
    'active_RMS_cartesian_kcal_mol_A': .03,
    'bond_change_A': 1e-9, 'angle_cosine_change': 1e-10,
    'protected_distance_change_A': 1e-9, 'fixed_atom_change_A': 1e-12,
    'displacement_roundoff_A': 1e-12}
NAMES = ['metal_x', 'metal_y', 'metal_z', 'A7:sidechain:CA-CB',
         'A7:sidechain:CB-CG', 'A3:sidechain:CA-CB', 'A11:sidechain:CA-CB']


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
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def relative(value, origin):
    mace = math.fsum([value['MACE_eV'], -origin['MACE_eV']]) * EV_TO_KCAL
    water = math.fsum([value['GFN2_water_hartree'], -origin['GFN2_water_hartree']]) * HA_TO_KCAL
    vacuum = math.fsum([value['GFN2_vacuum_hartree'], -origin['GFN2_vacuum_hartree']]) * HA_TO_KCAL
    return {'MACE_kcal_mol': mace, 'GFN2_water_kcal_mol': water, 'GFN2_vacuum_kcal_mol': vacuum,
            'solvent_kcal_mol': math.fsum([water, -vacuum]),
            'composite_kcal_mol': math.fsum([mace, water, -vacuum])}


class GeometryGuard:
    """Original peptide guards plus all-atom moving-metal contacts and bounds."""
    def __init__(self, source, active, Geometry, settings):
        self.source = source; self.x = source['x']; self.atoms = source['atoms']; self.active = active
        self.kin = source['kin']; self.settings = settings
        self.bonds = [b['native_particle_indices'] for b in source['bonds']]
        self.moving = set(i for k in active for i in self.kin.modes[k]['moving_indices'])
        require(255 in self.moving, 'Active target translations missing')
        self.fixed = sorted(set(range(256)) - self.moving)
        adapted = [dict(source_id=a['id'], residue_id=a['id'].rsplit('/', 1)[0],
            name=a['id'].rsplit('/', 1)[1], element=a['element'], resname=a['resname']) for a in self.atoms[:255]]
        self.peptide = Geometry(adapted, self.bonds, self.x[:255], self.x[:255],
            sorted(self.moving - {255}), [{'element': METAL, 'xyz_A': self.x[255].tolist()}])
        graph = [set() for _ in self.atoms]; self.excluded = set()
        for a, b in self.bonds:
            graph[a].add(b); graph[b].add(a)
        for i, neighbors in enumerate(graph):
            near = set(neighbors)
            for j in neighbors:
                near.update(graph[j])
            self.excluded.update(tuple(sorted((i, j))) for j in near if i != j)
        self.angles = np.asarray([(a, b, c) for b, ns in enumerate(graph)
                                  for a in sorted(ns) for c in sorted(ns) if a < c])
        self.original_cosines = self.cosines(self.x)
        self.groups = [g['indices'] for g in source['modes']['protected_groups']]

    def cosines(self, positions):
        a, b, c = self.angles.T; u, v = positions[a] - positions[b], positions[c] - positions[b]
        return np.einsum('ij,ij->i', u, v) / np.linalg.norm(u, axis=1) / np.linalg.norm(v, axis=1)

    def check(self, q, positions, previous):
        from scipy.spatial import cKDTree
        require(positions.shape == (256, 3) and np.isfinite(positions).all(), 'Invalid candidate coordinates')
        self.peptide.other = positions[255:].copy()  # all_source remains the original immutable baseline
        inherited = self.peptide.check(positions[:255]); overlaps = []
        for a, b in sorted(cKDTree(positions).query_pairs(self.settings['new_overlap_heavy_A'])):
            if (a not in self.moving and b not in self.moving) or (a, b) in self.excluded:
                continue
            cutoff = self.settings['new_overlap_H_A'] if 'H' in (self.atoms[a]['element'], self.atoms[b]['element']) else self.settings['new_overlap_heavy_A']
            distance = float(np.linalg.norm(positions[a] - positions[b]))
            original = float(np.linalg.norm(self.x[a] - self.x[b]))
            if distance < min(cutoff, original) - self.settings['overlap_roundoff_A']:
                overlaps.append({'indices': [a, b], 'source_ids': [self.atoms[a]['id'], self.atoms[b]['id']],
                                 'source_A': original, 'candidate_A': distance, 'cutoff_A': cutoff})
        extent = float(np.linalg.norm(positions - self.x, axis=1).max())
        step = float(np.linalg.norm(positions - previous, axis=1).max())
        target = float(np.linalg.norm(positions[255] - self.x[255]))
        torsion = max(abs(float(q[k])) for k in self.active if self.kin.modes[k]['unit'] == 'radian')
        fixed = float(np.max(np.abs(positions[self.fixed] - self.x[self.fixed]), initial=0))
        angle = float(np.max(np.abs(self.cosines(positions) - self.original_cosines)))
        protected = max(float(np.max(np.abs(np.linalg.norm(positions[g][:, None] - positions[g][None, :], axis=2) -
                                            np.linalg.norm(self.x[g][:, None] - self.x[g][None, :], axis=2)))) for g in self.groups)
        eps = SETTINGS['displacement_roundoff_A']
        passed = bool(inherited['pass'] and not overlaps and extent <= .20 + eps and target <= .20 + eps
                      and step <= .05 + eps and torsion <= .20 + 1e-12 and fixed <= 1e-12
                      and inherited['source_bond_max_error_A'] <= 1e-9 and angle <= 1e-10 and protected <= 1e-9)
        return {'pass': passed, 'peptide_guard': inherited, 'all256_new_worsened_severe_contacts': overlaps,
            'maximum_source_displacement_A': extent, 'target_displacement_A': target,
            'maximum_trial_displacement_A': step, 'maximum_torsion_radian': torsion,
            'fixed_atom_max_change_A': fixed, 'bond_angle_cosine_max_change': angle,
            'protected_group_distance_max_change_A': protected,
            'boundary_flags': {'source_displacement': extent >= .20 - 1e-5,
                'target_displacement': target >= .20 - 1e-5, 'torsion': torsion >= .20 - 1e-5}}


class CompositeEvaluator:
    """Current native MACE helper and the already admitted standalone gradient path."""
    def __init__(self, world, source, out):
        self.world, self.source, self.out = world, source, out
        self.native = None; self.cache = {}; self.records = []; self.requests = []
        self.lock = threading.Lock()
        self.cost = {'new_composite_evaluations': 0, 'MACE_calls_attempted': 0,
                     'GFN2_calls_attempted': 0, 'MACE_call_seconds': 0., 'GFN2_call_seconds': 0.,
                     'model_load_seconds': 0., 'cache_hits': 0}
        self.base = source['prepared']['endpoints']['whole'][METAL]
        self.symbols = [METAL if a['element'] == 'M' else a['element'] for a in source['atoms']]
        self.recipe = {'protocol': PROTOCOL, 'preparation': world['modes']['preparation'],
            'model': world['mcfg']['model'], 'MACE_versions': world['mcfg']['expected_versions'],
            'GFN2_runtime': world['scfg']['runtime'], 'GFN2_settings': world['solvent_module'].EXPECTED_SETTINGS,
            'GFN2_tolerances': world['scfg']['tolerances'], 'physical_charge': -2,
            'physical_multiplicity': 1, 'GFN2_effective_multiplicity': 1,
            'physical_ids': source['modes']['physical_source_ids'], 'symbols': self.symbols}
        for name in ('executable', 'parameters'):
            verify(world['scfg']['runtime'][name])
        for record in world['scfg']['runtime']['runtime_libraries']:
            verify(record)
        self.solvent = world['solvent_module']
        self.backend = module('lbt3_accommodation_backend', world['scfg']['backend_helper'])
        self.gradients = module('lbt3_accommodation_gradient_parser', world['scfg']['gradient_helper'])
        self.common = module('lbt3_accommodation_common', world['scfg']['common_helper'])
        require(self.common.HA_TO_KCAL == HA_TO_KCAL and self.common.BOHR_TO_A == .529177210903, 'Unit convention differs')
        self.raw = world['raw']
        self.origin = self.save_origin()

    def key(self, positions):
        return digest({'recipe': self.recipe, 'positions_A': np.asarray(positions).tolist()})

    def save_origin(self):
        directory = self.out / 'origin'; directory.mkdir()
        np.save(directory / 'POSITIONS.npy', self.source['x'], allow_pickle=False)
        np.save(directory / 'GRADIENT.npy', -self.source['forces']['composite'][METAL], allow_pickle=False)
        raw = self.world['force_by'][self.base['task_id']]
        components = {'MACE_eV': raw['energy_eV'],
            'GFN2_vacuum_hartree': self.raw[(self.base['task_id'], 'vacuum')]['energy_hartree'],
            'GFN2_water_hartree': self.raw[(self.base['task_id'], 'ALPB(Water)')]['energy_hartree']}
        point = {'status': 'complete', 'cache_key': self.key(self.source['x']), 'recipe': self.recipe,
            'source_id': SOURCE, 'metal': METAL, 'components': components,
            'positions': pin(directory / 'POSITIONS.npy'), 'gradient_kcal_mol_A': pin(directory / 'GRADIENT.npy'),
            'xyz': self.base['xyz'], 'physical_source': self.base['physical_source'],
            'reused_origin': {'MACE': self.world['profile']['mace_forces'],
                              'solvent': self.world['profile']['solvent_result'], 'origin_task_id': self.base['task_id']},
            'new_MACE_calls': 0, 'new_GFN2_calls': 0}
        write(directory / 'POINT.json', point); record = pin(directory / 'POINT.json')
        self.cache[point['cache_key']] = record; self.records.append(record)
        return point, record

    def low(self, parent, medium, physical_pin):
        old = self.raw[(self.base['task_id'], medium)]['task']
        directory = parent / ('water' if medium == 'ALPB(Water)' else 'vacuum'); directory.mkdir()
        (directory / 'core.xyz').write_bytes((parent / 'whole.xyz').read_bytes())
        task = dict(old, task_id=parent.name + '_' + METAL + '_' + directory.name,
            directory=str(directory), xyz=pin(directory / 'core.xyz'), source_xyz=pin(parent / 'whole.xyz'),
            physical_source=physical_pin, origin_task_id=self.base['task_id'])
        task.pop('task_record', None); write(directory / 'TASK.json', task)
        task['task_record'] = pin(directory / 'TASK.json')
        row = {'status': 'unavailable', 'task': task, 'call_attempted': False}
        tick = time.monotonic()
        try:
            env = os.environ.copy()
            for key, value in self.world['manifest']['effective_environment'].items():
                if value is None:
                    env.pop(key, None)
                else:
                    env[key] = value
            require(env.get('OMP_NUM_THREADS') == '1', 'Qualified standalone thread budget differs')
            with self.lock:
                self.cost['GFN2_calls_attempted'] += 1
            row['call_attempted'] = True
            execution = self.backend.call(task['command'], directory, env)
            row['execution'] = pin(directory / 'EXECUTION.json')
            verify(task['xyz'])
            row.update(self.solvent.admit(directory, task, execution, self.backend, self.gradients, self.common, np), status='accepted')
        except Exception as exc:
            row.update(failure=str(exc), traceback=traceback.format_exc())
        row['wall_seconds'] = time.monotonic() - tick
        with self.lock:
            self.cost['GFN2_call_seconds'] += row['wall_seconds']
        row['artifacts'] = {p.name: pin(p) for p in directory.iterdir() if p.is_file()}
        write(directory / 'RESULT.json', row)
        return row

    def evaluate(self, positions):
        key = self.key(positions)
        if key in self.cache:
            self.cost['cache_hits'] += 1; record = self.cache[key]; point = read(record)
            require(point['status'] == 'complete', 'Previously failed composite evaluation cannot be retried')
            self.requests.append({'point': record, 'reused': True})
            return point, record
        require(self.cost['new_composite_evaluations'] < 24, 'Composite call budget exceeded')
        self.cost['new_composite_evaluations'] += 1
        directory = self.out / ('evaluation_%03d' % self.cost['new_composite_evaluations']); directory.mkdir()
        tick = time.monotonic()
        point = {'status': 'unavailable', 'source_id': SOURCE, 'metal': METAL, 'cache_key': key,
                 'recipe': self.recipe, 'new_MACE_calls': 0, 'new_GFN2_calls': 0}
        np.save(directory / 'POSITIONS.npy', positions, allow_pickle=False)
        point['positions'] = pin(directory / 'POSITIONS.npy')
        physical = [dict(a, xyz_A=positions[i].tolist()) for i, a in enumerate(self.source['atoms'])]
        write(directory / 'PHYSICAL.json', physical); point['physical_source'] = pin(directory / 'PHYSICAL.json')
        with (directory / 'whole.xyz').open('x') as handle:
            handle.write('256\nLBT3 whole own-metal composite response\n')
            for symbol, p in zip(self.symbols, positions):
                handle.write(symbol + ' ' + ' '.join(format(v, '.17g') for v in p) + '\n')
        point['xyz'] = pin(directory / 'whole.xyz')
        try:
            if self.native is None:
                load_start = time.monotonic()
                self.native = self.world['fd_module'].NativeMACE(self.world['mcfg'])
                self.cost['model_load_seconds'] = time.monotonic() - load_start
            # Native evaluation is complete before two single-thread CPU workers start.
            start_calls = self.native.calls; mace_tick = time.monotonic()
            try:
                native = self.native.evaluate(self.base, self.symbols, positions, compute_force=True)
            finally:
                attempted = self.native.calls - start_calls
                self.cost['MACE_calls_attempted'] += attempted; point['new_MACE_calls'] = attempted
                self.cost['MACE_call_seconds'] += time.monotonic() - mace_tick
            require(native['force_convention'] == 'negative derivative of the same native vacuum energy', 'Native force convention differs')
            write(directory / 'MACE.json', native); point['MACE'] = pin(directory / 'MACE.json')
            before = self.cost['GFN2_calls_attempted']
            with ThreadPoolExecutor(max_workers=2) as pool:
                low = list(pool.map(lambda medium: self.low(directory, medium, point['physical_source']), ('vacuum', 'ALPB(Water)')))
            point['new_GFN2_calls'] = self.cost['GFN2_calls_attempted'] - before
            point['GFN2'] = {row['task']['medium']: pin(Path(row['task']['directory']) / 'RESULT.json') for row in low}
            require(all(row['status'] == 'accepted' for row in low), 'Required standalone energy/gradient unavailable')
            vacuum, water = low
            gradient = -np.asarray(native['forces_eV_A']) * EV_TO_KCAL + np.asarray(water['gradient']['kcal_mol_A']) - np.asarray(vacuum['gradient']['kcal_mol_A'])
            require(gradient.shape == (256, 3) and np.isfinite(gradient).all(), 'Composite gradient invalid')
            np.save(directory / 'GRADIENT.npy', gradient, allow_pickle=False)
            point.update(status='complete', gradient_kcal_mol_A=pin(directory / 'GRADIENT.npy'),
                components={'MACE_eV': native['energy_eV'], 'GFN2_vacuum_hartree': vacuum['energy_hartree'],
                            'GFN2_water_hartree': water['energy_hartree']})
        except Exception as exc:
            point.update(failure=str(exc), traceback=traceback.format_exc())
        point['elapsed_seconds'] = time.monotonic() - tick
        write(directory / 'POINT.json', point); record = pin(directory / 'POINT.json')
        self.cache[key] = record; self.records.append(record); self.requests.append({'point': record, 'reused': False})
        require(point['status'] == 'complete', 'Composite backend/state/gradient failure: ' + point.get('failure', 'unknown'))
        return point, record


def projected_diagnostic(jac, gradient, active, modes):
    full = jac.reshape(64, -1).T; scales = np.linalg.norm(full, axis=0)
    require(np.all(scales > 1e-12), 'Zero current physical mode')
    normed = full / scales; matrix = normed[:, active]; gram = matrix.T @ matrix
    eigenvalues, vectors = np.linalg.eigh(gram)
    keep = eigenvalues > 1e-10 * eigenvalues[-1]
    require(int(keep.sum()) == 7, 'Active metric rank changed')
    inverse = (vectors[:, keep] / eigenvalues[keep]) @ vectors[:, keep].T
    flat = gradient.reshape(-1); projected = (matrix @ inverse @ matrix.T @ flat).reshape(256, 3)
    residual = gradient - projected; loads = -normed.T @ flat
    record = {'active_span_max_atom_kcal_mol_A': float(np.linalg.norm(projected, axis=1).max()),
        'active_span_RMS_cartesian_kcal_mol_A': float(np.sqrt(np.mean(projected ** 2))),
        'active_span_L2_kcal_mol_A': float(np.linalg.norm(projected)),
        'full_gradient_L2_kcal_mol_A': float(np.linalg.norm(gradient)),
        'outside_active_span_L2_kcal_mol_A': float(np.linalg.norm(residual)),
        'outside_active_span_max_atom_kcal_mol_A': float(np.linalg.norm(residual, axis=1).max()),
        'normalized_Gram_eigenvalues': eigenvalues.tolist(), 'rank': 7,
        'mode_force_loads_kcal_mol_A': [{'mode_index': k, 'name': mode['name'], 'family': mode['family'],
            'active': k in active, 'force_load': float(loads[k])} for k, mode in enumerate(modes)],
        'stationary_minimum_claimed': False}
    record['low_active_span_force'] = bool(record['active_span_max_atom_kcal_mol_A'] <= .1
                                          and record['active_span_RMS_cartesian_kcal_mol_A'] <= .03)
    return record


def search(source, active, evaluator, out, Geometry, geometry_settings, lbfgs_direction):
    kin = source['kin']; x = source['x']; n = 7
    q0 = np.zeros(64); replay, _, jac0, _ = kin.evaluate(q0)
    matrix = jac0[active].reshape(n, -1).T; scales = np.linalg.norm(matrix, axis=0)
    normalized = matrix / scales; eig, vectors = np.linalg.eigh(normalized.T @ normalized)
    require(int(np.sum(eig > 1e-10 * eig[-1])) == 7, 'Source active physical metric is rank deficient')
    inverse_sqrt = (vectors / np.sqrt(eig)) @ vectors.T
    transform = inverse_sqrt / scales[:, None]
    require(np.allclose((matrix @ transform).T @ (matrix @ transform), np.eye(7), atol=1e-10, rtol=0), 'Whitening replay failed')
    guard = GeometryGuard(source, active, Geometry, geometry_settings)
    origin_guard = guard.check(q0, x, x); require(origin_guard['pass'], 'Actual origin fails geometry guards')
    result = {'source_id': SOURCE, 'metal': METAL, 'status': 'unavailable', 'candidate_admitted': False,
        'active_indices': active, 'active_names': NAMES, 'coordinate_scaling_T': transform.tolist(),
        'source_column_norms': scales.tolist(), 'source_normalized_Gram_eigenvalues': eig.tolist(),
        'q0_kinematic_replay_max_error_A': float(np.max(np.abs(replay - x))),
        'authoritative_origin': source['prepared']['endpoints']['whole'][METAL]['xyz'],
        'origin_coordinate_policy': 'exact physical source and saved energy/force; q0 replay roundoff never triggers new evaluation',
        'origin_guard': origin_guard, 'trace': [], 'accepted_steps': 0, 'stationary_minimum_claimed': False}
    require(result['q0_kinematic_replay_max_error_A'] <= 1e-12, 'Origin kinematic replay differs')
    u = np.zeros(7); q = q0; positions = x.copy(); jac = jac0; history = []
    point, point_pin = evaluator.origin; origin = point
    gradient = np.load(verify(point['gradient_kcal_mol_A']), allow_pickle=False)
    gu = transform.T @ np.einsum('kij,ij->k', jac[active], gradient)
    current_guard = origin_guard
    try:
        while True:
            diagnostic = projected_diagnostic(jac, gradient, active, kin.modes)
            row = {'accepted_step': result['accepted_steps'], 'u_A': u.tolist(), 'full_q': q.tolist(),
                   'point': point_pin, 'work_from_origin': relative(point['components'], origin['components']),
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
            for multiplier in SETTINGS['trial_multipliers']:
                trial_u = u + multiplier * direction; trial_q = np.zeros(64); trial_q[active] = transform @ trial_u
                trial_positions, _, trial_jac, _ = kin.evaluate(trial_q)
                check = guard.check(trial_q, trial_positions, positions)
                trial = {'multiplier': multiplier, 'u_A': trial_u.tolist(), 'full_q': trial_q.tolist(), 'geometry': check,
                         'molecular_evaluation_requested': False, 'accepted': False}
                row['trials'].append(trial)
                if not check['pass']:
                    trial['reason'] = 'geometry_or_bound'; continue
                if evaluator.cost['new_composite_evaluations'] >= 24 and evaluator.key(trial_positions) not in evaluator.cache:
                    trial['reason'] = 'evaluation_limit'; break
                trial['molecular_evaluation_requested'] = True
                next_point, next_pin = evaluator.evaluate(trial_positions)
                change = relative(next_point['components'], point['components'])
                trial.update(point=next_pin, work_from_current=change, armijo_rhs_kcal_mol=1e-4 * multiplier * slope)
                if change['composite_kcal_mol'] <= trial['armijo_rhs_kcal_mol']:
                    trial['accepted'] = True
                    accepted = (trial_u, trial_q, trial_positions, trial_jac, next_point, next_pin, check)
                    break
                trial['reason'] = 'own_composite_Armijo'
            # Save completed rejected trials even if the finite line search stops.
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
        work = relative(point['components'], origin['components'])
        require(work['composite_kcal_mol'] <= 1e-8, 'Last accepted own energy increased')
        result.update(status='finite_proposal_available', candidate_admitted=True, termination=termination,
            candidate=point_pin, full_q=q.tolist(), u_A=u.tolist(), final_geometry=final_check,
            work_from_origin=work, final_force_diagnostic=diagnostic,
            preserved_boundary_flags=final_check['boundary_flags'])
    except Exception as exc:
        result.update(status='unavailable', candidate_admitted=False, termination='backend_or_numerical_failure',
                      failure=str(exc), traceback=traceback.format_exc(), last_archived_valid_point=point_pin)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('config', 'output'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--repo-root', type=Path, default=REPO)
    parser.add_argument('--stage', choices=['scout'], default='scout')
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    cpus = int(os.environ['SLURM_CPUS_PER_TASK']); require(cpus >= 2, 'Two single-thread solvent workers require two CPUs')
    require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'), 'Allocated GPU required')
    cfg = json.loads(args.config.read_text()); verify(cfg['plan'])
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Adapter differs from frozen config')
    profile = read(cfg['profile_result']); fd_result = read(cfg['fd_result'])
    require(profile['status'] == 'complete' and fd_result['status'] == 'passed' and fd_result['composite_passed_all_states']
            and fd_result['design']['profile_result'] == cfg['profile_result'], 'Matching complete derivative gate required')
    require(profile['active_candidate_basis']['mode_names_selection_order'] == NAMES, 'Frozen seven-coordinate basis differs')
    verify(fd_result['config']); fd_module = module('lbt3_accommodation_native', fd_result['implementation'])
    profile_module = module('lbt3_accommodation_profile', profile['implementation'])
    modes = read(profile['modes_result']); solvent = read(profile['solvent_result'])
    mforce = read(profile['mace_forces']); mcfg = read(mforce['config']); scfg = read(solvent['config'])
    require(modes['status'] == solvent['status'] == mforce['status'] == 'complete'
            and modes['preparation'] == solvent['preparation'] == mcfg['preparation'], 'Origin source lineage differs')
    require(mforce['model'] == mcfg['model'] and mforce['versions'] == mcfg['expected_versions'], 'Origin native model/runtime differs')
    solvent_module = module('lbt3_accommodation_solvent', solvent['implementation'])
    require(scfg['tolerances'] == solvent_module.TOLERANCES, 'Solvent admission settings differ')
    for name, record in profile['helpers'].items():
        verify(record); require(pin(args.repo_root / 'scripts' / name)['sha256'] == record['sha256'], 'Helper changed: ' + name)
    sys.path.insert(0, str(args.repo_root.resolve() / 'scripts'))
    global np
    import numpy as np
    from mace_site_kinematics import Kinematics
    from scaffold_accommodation import Geometry, SETTINGS as GEOMETRY_SETTINGS, lbfgs_direction
    prep = read(modes['preparation']); mode_rows = {r['source_id']: r for r in modes['rows']}
    prep_rows = {r['source_id']: r for r in prep['rows']}
    corrections = {(r['source_id'], r['metal']): r for r in solvent['corrections']}
    raw = {(r['task']['origin_task_id'], r['task']['medium']): r for r in solvent['rows']}
    force_by = {r['task_id']: r for r in mforce['rows']}
    require(len(raw) == 8 and len(force_by) == len(corrections) == 4
            and all(r['status'] == 'accepted' for r in raw.values()), 'Whole origin matrix incomplete')
    source = profile_module.source_profile(SOURCE, mode_rows[SOURCE], prep_rows[SOURCE], corrections,
        force_by, raw, solvent['gradient_units'], np, Kinematics)
    active = profile['active_candidate_basis']['mode_indices_selection_order']
    require([source['kin'].modes[k]['name'] for k in active] == NAMES, 'Active indices differ from physical maps')
    world = dict(profile=profile, modes=modes, mcfg=mcfg, scfg=scfg, raw=raw, force_by=force_by,
                 solvent_module=solvent_module, fd_module=fd_module, manifest=read(solvent['manifest']))
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); start = time.monotonic()
    result = {'protocol_id': PROTOCOL, 'stage': 'scout', 'status': 'unavailable', 'settings': SETTINGS,
        'config': pin(args.config), 'implementation': pin(__file__), 'plan': cfg['plan'],
        'profile_result': cfg['profile_result'], 'fd_result': cfg['fd_result'],
        'native_helper': fd_result['implementation'], 'solvent_helper': solvent['implementation'],
        'profile_helper': profile['implementation'], 'helpers': profile['helpers'],
        'source_id': SOURCE, 'metal': METAL, 'whole_peptide_atoms': 256, 'charge': -2,
        'physical_multiplicity': 1, 'GFN2_effective_multiplicity': 1,
        'added_forcefield_energy': False, 'stationary_minimum_claimed': False,
        'pool_comparison_performed': False, 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname()}
    write(out / 'START.json', result)
    evaluator = None
    try:
        evaluator = CompositeEvaluator(world, source, out)
        row = search(source, active, evaluator, out, Geometry, GEOMETRY_SETTINGS, lbfgs_direction)
        write(out / 'SEARCH.json', row); result['search'] = pin(out / 'SEARCH.json')
        result['status'] = 'complete' if row['candidate_admitted'] else 'unavailable'
        result['continuation_execution_and_physical_admission'] = row['candidate_admitted']
        result['continuation_rule'] = 'root inspects scout before unchanged other three; no sign/improvement criterion'
    except Exception as exc:
        result.update(failure=str(exc), traceback=traceback.format_exc())
    if evaluator is not None:
        result.update(actual_cost=evaluator.cost, points=evaluator.records, requests=evaluator.requests)
        if evaluator.native is not None:
            result['native_runtime'] = {'versions': evaluator.native.versions, 'device': evaluator.native.device,
                'actual_model_forwards': evaluator.native.calls, 'torch_threads': evaluator.native.torch.get_num_threads(),
                'peak_cuda_allocated_bytes': evaluator.native.torch.cuda.max_memory_allocated()}
    result['elapsed_seconds'] = time.monotonic() - start
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds'] * cpus
    result['peak_host_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['allocation'] = {k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_CPU', 'CUDA_VISIBLE_DEVICES')}
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}), flush=True)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
