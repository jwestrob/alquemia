#!/usr/bin/env python3
"""Whole256 LBT Cartesian response with the unchanged admitted composite."""
import argparse
import copy
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
import traceback

PROTOCOL = 'lbt_full_cartesian_matched_energy_v1'
STATES = (('7CCO', 'La'), ('7CCO', 'Dy'), ('7CCN', 'La'), ('7CCN', 'Dy'))
FD_STEPS = (1e-4, 5e-5)
FD_ABSOLUTE = .005
FD_RELATIVE = 1e-5
LIMITS = {'scout': {'accepted_steps': 12, 'evaluations': 40},
          'run': {'accepted_steps': 80, 'evaluations': 200}}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    p = Path(record['path'])
    require(pin(p)['sha256'] == record['sha256'], 'Changed pin: ' + str(p))
    return p


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def load_inputs(cfg, out):
    """Saved replay only. Returns base/world/terminals/guards for root's pool."""
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required, including preflight')
    verify(cfg['implementation']); verify(cfg['plan']); verify(cfg['spec'])
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Cartesian adapter differs')
    core = module('lbt_cartesian_core', cfg['optimizer'])
    parent = read(cfg['parent_pool_result']); pcfg = read(parent['config'])
    verify(parent['implementation']); verify(pcfg['implementation'])
    require(parent['status'] == 'complete' and parent['slurm_job_id'] == '12767'
            and parent['implementation']['sha256'] == pcfg['implementation']['sha256'] ==
            '6b4944ac75918816fb303ba2e7ce1a504710303240e448ce1ca522f1c2dc2a3a'
            and parent['declared_cells'] == 44 and parent['declared_geometries'] == 22,
            'Exact completed parent pool12767 required')
    pool = module('lbt_cartesian_parent_pool', pcfg['implementation'])
    admission = out/'parent_admission'; admission.mkdir()
    replay=module('lbt_cartesian_replay_v2',cfg['replay_helper'])
    base, world, prior, candidates, cells, slots, parents = replay.pool_inputs(pool,pcfg,admission,cfg['recovery_evidence'])
    admitted, _ = pool.assemble(base, world, candidates, cells, slots)
    require(len(admitted) == 4 and all(r['status'] == 'admitted' for r in admitted), 'Four prior own terminals required')
    require(parents[0]['slurm_job_id'] == '12762' and parents[1]['slurm_job_id'] == '12765', 'Exact parent searches required')
    global np
    import numpy as np
    geometry = importlib.import_module('scaffold_accommodation')
    terminals = {}; guards = {}; preflight = []
    actual_pool_cells = {(r['candidate_id'], m): c for r in parent['paired_matrix'] for m, c in r['cells'].items()}
    for row in admitted:
        sid, metal = row['source_id'], row['metal']; state = (sid, metal)
        require(state in STATES and state not in terminals, 'Duplicate/unknown state')
        saved = next(r for r in parent['terminal_slots'] if (r['source_id'], r['metal']) == state)
        require(all(saved[k] == row[k] for k in ('status', 'candidate_id', 'point', 'search', 'alias')),
                'Saved terminal identity differs')
        source = world['sources'][sid]; atoms = source['atoms']; ids = [a['id'] for a in atoms]
        point, x, native = base.admit_point(world, source, metal, row['point'], True)
        require(actual_pool_cells[(row['candidate_id'], metal)]['point'] == row['point']
                and actual_pool_cells[(row['candidate_id'], metal)]['native_energies'] == native,
                'Paired pool own cell differs')
        require(len(atoms) == 256 and ids[-1] == 'metal' and source['prepared']['waters'] ==
                source['prepared']['spectators'] == source['prepared']['caps'] == [], 'Whole LBT inventory differs')
        bond_rows = source['bonds']; pairs = [b['native_particle_indices'] for b in bond_rows]
        for b, pair in zip(bond_rows, pairs):
            require([ids[i] for i in pair] == [b['atom_a_id'], b['atom_b_id']], 'Native bond physical mapping differs')
        eq = [b['forcefield_equilibrium_A'] for b in bond_rows]
        donor_map = read(source['prepared']['endpoints']['matched_core'][metal]['source_cap_map'])
        donor_ids = [f"{c['chain']}/{c['resnum']}//{c['atom']}" for c in donor_map['contacts']]
        expected = ['A/3//OD1', 'A/3//OD2', 'A/5//OD1', 'A/7//OD1', 'A/7//OD2',
                    'A/9//O', 'A/11//OE1', 'A/11//OE2', 'A/14//OE1', 'A/14//OE2']
        require(len(donor_ids) == 10 and set(donor_ids) == set(expected), 'Declared site oxygen inventory differs')
        target = {'index': 255, 'label': 'LBT_target', 'site_donor_indices': [ids.index(a) for a in expected],
                  'coordinating_water_indices': []}
        guard = core.PhysicalGuard(atoms, pairs, eq, source['x'], geometry, list(range(255)), [target])
        original_check = guard.check(source['x']); check = guard.check(x)
        require(original_check['pass'] and original_check['state_compatible'] and check['pass']
                and check['state_compatible'], 'Original/terminal physical admission unavailable')
        cache = dict(world['starts'][state]['cache']); records = list(world['starts'][state]['records'])
        owner = next(r for r in slots if (r['source_id'], r['metal']) == state)
        for record in owner['points']:
            p = read(record)
            if 'gradient_kcal_mol_A' in p:
                p, _, _ = base.admit_point(world, source, metal, record, True)
                cache[p['cache_key']] = record
                if record not in records:
                    records.append(record)
        cache[point['cache_key']] = row['point']
        key = base.key_for(world, source, metal, x)
        changed = x.copy(); changed[0, 0] = np.nextafter(changed[0, 0], math.inf)
        recipe = base.recipe(world, source, metal); altered = copy.deepcopy(recipe)
        altered['protocol_id'] += ':not_an_executable_recipe'
        protocol_key = base.digest({'recipe': altered, 'positions_A': x.tolist()})
        other = 'Dy' if metal == 'La' else 'La'
        tests = {'coordinate_changes_key': base.key_for(world, source, metal, changed) != key,
                 'state_changes_key': base.key_for(world, source, other, x) != key,
                 'protocol_changes_key': protocol_key != key,
                 'exact_coordinate_replay_key': key == point['cache_key']}
        require(all(tests.values()), 'Coordinate/state/protocol cache invalidation failed')
        terminals[state] = {'source': source, 'point': point, 'point_pin': row['point'], 'search': row['search'],
                            'x': x, 'g': np.load(verify(point['gradient_kcal_mol_A']), allow_pickle=False),
                            'cache': cache, 'records': records, 'owner': owner, 'candidate_id': row['candidate_id']}
        guards[state] = guard
        preflight.append({'source_id': sid, 'metal': metal, 'pass': True, 'terminal_point': row['point'],
            'previous_search': row['search'], 'physical_charge': -2, 'physical_multiplicity': 1 if metal == 'La' else 6,
            'original_source': source['prepared']['physical_atoms'], 'original_guard': original_check,
            'terminal_guard': check, 'cache_checks': tests, 'full_gradient_replayed': True,
            'bond_count': len(pairs), 'all_atoms_movable': 256,
            'actual_full_force': core.force_metrics(terminals[state]['g']),
            'old_coordinate_constraints_used_for_new_trials': False})
    require(set(terminals) == set(STATES), 'Four-state preflight denominator differs')
    return {'base': base, 'world': world, 'terminals': terminals, 'guards': guards,
            'parent_pool': parent, 'parent_config': pcfg, 'core': core,
            'preflight': {'status': 'passed', 'declared_states': 4, 'new_molecular_calls': 0,
                          'model_initialized': False, 'rows': preflight},
            'geometry_helper': pin(geometry.__file__),
            'archived_metric_replay':world['preflight']['rows']}


class Evaluator:
    """Thin full-Cartesian view of the immutable matched-energy evaluator."""
    def __init__(self, loaded, state, out, native=None, additional_points=()):
        self.loaded = loaded; self.state = state; self.base = loaded['base']; self.world = loaded['world']
        self.start = loaded['terminals'][state]; self.source = self.start['source']; self.metal = state[1]
        out.mkdir()
        seed = {'cache': dict(self.start['cache']), 'records': list(self.start['records'])}
        for record in additional_points:
            p, _, _ = self.base.admit_point(self.world, self.source, self.metal, record, True)
            seed['cache'][p['cache_key']] = record
            if record not in seed['records']:
                seed['records'].append(record)
        self.engine = self.base.MatchedEvaluator(self.world, self.source, self.metal, out, seed, native=native)
        self.pin = pin

    @property
    def new_evaluations(self):
        return self.engine.cost['new_composite_evaluations']

    @property
    def evaluation_requests(self):
        return self.engine.requests

    def read_point(self, record):
        p, x, native = self.base.admit_point(self.world, self.source, self.metal, record, True)
        g = np.load(verify(p['gradient_kcal_mol_A']), allow_pickle=False)
        parts = {'MACE_kcal_mol': native['MACE_eV']*self.base.EV,
                 'fixed_ALPB_kcal_mol': native['ALPB_Eh']*self.base.HA}
        return {'total_energy_kcal_mol': math.fsum(parts.values()), 'positions_A': x,
                'gradient_kcal_mol_A': g, 'component_energies_kcal_mol': parts,
                'native_components': p['components'], 'provenance': record, 'physical_state': p['recipe']}

    def evaluate(self, positions):
        _, record, _ = self.engine.evaluate(positions, compute_force=True)
        return self.read_point(record)

    def difference(self, new, old):
        return self.base.relative(new['native_components'], old['native_components'])

    def gradients(self, record):
        point = read(record); ml = read(point['MACE']); sol = read(point['ALPB'])
        mace = -np.asarray(ml['forces_eV_A'])*self.base.EV
        alpb = np.asarray([[float(v) for v in line.split()] for line in
                           verify(sol['artifacts']['gradient.dat']).read_text().splitlines()])*(self.base.HA/self.base.BOHR)
        require(np.array_equal(mace+alpb, self.read_point(record)['gradient_kcal_mol_A']), 'Component gradient closure differs')
        return {'MACE_kcal_mol': mace, 'fixed_ALPB_kcal_mol': alpb, 'composite_kcal_mol': mace+alpb}


def directions(loaded, state):
    start = loaded['terminals'][state]; x = start['x']; ids = [a['id'] for a in start['source']['atoms']]
    raw = np.sin((np.arange(768).reshape(256, 3)+1)*np.sqrt(2.))
    coupled, _ = loaded['core'].remove_rigid(x, raw)
    i, j = [ids.index(a) for a in ('A/8//C', 'A/9//N')]
    require(any(set(b['native_particle_indices']) == {i, j} for b in start['source']['bonds']),
            'Declared finite-difference backbone bond absent')
    axis = x[j]-x[i]; axis /= np.linalg.norm(axis)
    stretch = np.zeros_like(x); stretch[i] = -axis; stretch[j] = axis
    return [('coupled_all256', coupled/np.linalg.norm(coupled, axis=1).max()),
            ('A8C_A9N_stretch', stretch/np.linalg.norm(stretch, axis=1).max())]


def qualify(cfg, loaded, out, result):
    native = None; rows = []; costs = []; stopped = False
    for state in STATES:
        directory = out/('_'.join(state)); directory.mkdir()
        evaluator = Evaluator(loaded, state, directory/'evaluations', native=native)
        original = evaluator.read_point(loaded['terminals'][state]['point_pin'])
        component_g = evaluator.gradients(original['provenance'])
        for name, direction in directions(loaded, state):
            row = {'source_id': state[0], 'metal': state[1], 'direction': name, 'status': 'unavailable',
                   'direction_vectors': direction.tolist(), 'origin': original['provenance'],
                   'cells': [{'step_A': sign*h, 'status': 'unavailable', 'reason': 'not attempted'}
                             for h in FD_STEPS for sign in (-1, 1)], 'checks': []}
            rows.append(row)
            if stopped:
                row['reason'] = 'not attempted after earlier unavailable qualification'; continue
            try:
                analytic = {k: float(np.sum(g*direction)) for k, g in component_g.items()}
                slopes = []
                for h in FD_STEPS:
                    pair = {}
                    for sign in (-1, 1):
                        x = original['positions_A']+sign*h*direction
                        gate = loaded['guards'][state].check(x)
                        cell = next(c for c in row['cells'] if c['step_A'] == sign*h)
                        cell.update(geometry=gate, reason='physical or execution admission pending')
                        require(gate['pass'] and gate['state_compatible'], 'Declared finite point fails physical admission')
                        point = evaluator.evaluate(x); pair[sign] = point
                        cell.update(status='complete', point=point['provenance']); cell.pop('reason')
                    difference = evaluator.difference(pair[1], pair[-1])
                    derivative = {k: difference[k]/(2*h) for k in analytic}; slopes.append(derivative)
                    for k, a in analytic.items():
                        tolerance = FD_ABSOLUTE+FD_RELATIVE*abs(a)
                        error = abs(derivative[k]-a)
                        row['checks'].append({'component': k, 'step_A': h, 'analytic': a, 'finite_difference': derivative[k],
                                              'absolute_error': error, 'tolerance': tolerance, 'pass': error <= tolerance})
                row['step_agreement'] = {k: {'absolute_difference': abs(slopes[0][k]-slopes[1][k]),
                    'tolerance': FD_ABSOLUTE+FD_RELATIVE*abs(analytic[k]),
                    'pass': abs(slopes[0][k]-slopes[1][k]) <= FD_ABSOLUTE+FD_RELATIVE*abs(analytic[k])} for k in analytic}
                row['status'] = 'passed' if all(c['pass'] for c in row['checks']) and all(c['pass'] for c in row['step_agreement'].values()) else 'failed_derivative'
                if row['status'] != 'passed':
                    stopped = True
            except Exception as exc:
                row.update(reason=repr(exc), traceback=traceback.format_exc()); stopped = True
        native = evaluator.engine.native
        costs.append(dict(source_id=state[0], metal=state[1], **evaluator.engine.cost))
    result.update(qualification_rows=rows, actual_costs=costs, declared_directions=8, declared_cells=32,
                  status='complete' if len(rows) == 8 and all(r['status'] == 'passed' for r in rows) else 'unavailable',
                  qualification_passed=len(rows) == 8 and all(r['status'] == 'passed' for r in rows),
                  derivative_steps_A=list(FD_STEPS), derivative_tolerance={'absolute': FD_ABSOLUTE, 'relative': FD_RELATIVE})
    require(sum(c['new_composite_evaluations'] for c in costs) <= 32, 'Qualification call ceiling exceeded')


def parent_gate(cfg, key, stage):
    r = read(cfg[key]); rcfg = read(r['config'])
    verify(r['implementation']); verify(rcfg['implementation'])
    require(r['status'] == 'complete' and r['stage'] == stage and r['protocol_id'] == PROTOCOL
            and r['implementation']['sha256'] == rcfg['implementation']['sha256'] == cfg['implementation']['sha256'],
            'Required Cartesian parent stage unavailable')
    for k in ('plan', 'spec', 'optimizer', 'parent_pool_result', 'replay_helper', 'recovery_evidence'):
        verify(rcfg[k]); verify(cfg[k])
        require(rcfg[k]['sha256'] == cfg[k]['sha256'], 'Cartesian parent recipe differs: ' + k)
    return r


def run_search(cfg, loaded, out, result, state, stage):
    qualification = parent_gate(cfg, 'qualification_result', 'qualify')
    require(qualification['qualification_passed'] and qualification['declared_directions'] == 8
            and qualification['declared_cells'] == 32
            and all(r['status'] == 'passed' for r in qualification['qualification_rows']), 'Complete finite derivative gate required')
    resume = None; extra = []; scout_row = None
    if stage == 'run':
        scout = parent_gate(cfg, 'scout_result', 'scout')
        require(scout['qualification_result'] == cfg['qualification_result'] and len(scout['rows']) == 1,
                'Scout qualification lineage differs')
        scout_row = scout['rows'][0]; sr = read(scout_row['search'])
        require((scout_row['source_id'], scout_row['metal']) == STATES[0] and scout_row['status'] == 'complete'
                and scout_row['physical_state_admitted'] and sr['physical_state_admitted']
                and sr['termination'] != 'execution_failure', 'Scout execution/physical gate failed')
        if state == STATES[0]:
            resume = read(sr['checkpoint']); extra = scout_row['points']
            require(resume['point'] == sr['terminal_point'] and resume['initial_point'] == sr['start_point']
                    == loaded['terminals'][state]['point_pin'] and resume['accepted_steps'] == sr['accepted_steps']
                    and resume['distinct_evaluations'] == sr['distinct_evaluations']
                    and resume['settings'] == loaded['core'].SETTINGS and resume['limits'] == LIMITS['scout'],
                    'Exact scout continuation checkpoint differs')
            require(0 <= resume['accepted_steps'] <= 12 and 0 <= resume['distinct_evaluations'] <= 40
                    and loaded['core'].SETTINGS['minimum_trust_A'] <= resume['trust_A'] <= loaded['core'].SETTINGS['maximum_trust_A'],
                    'Scout counters/trust invalid')
            for h in resume['history']:
                require(np.asarray(h['s']).shape == np.asarray(h['y']).shape == (768,)
                        and np.isfinite(h['s']).all() and np.isfinite(h['y']).all()
                        and math.isfinite(h['rho']) and h['rho'] > 0, 'Invalid scout L-BFGS history')
    d = out/('_'.join(state)); d.mkdir()
    evaluator = Evaluator(loaded, state, d/'evaluations', additional_points=extra)
    start = evaluator.read_point(loaded['terminals'][state]['point_pin'])
    search = loaded['core'].search(evaluator, loaded['guards'][state], start, d/'search', LIMITS[stage],
                                  loaded['world']['lbfgs'], resume=resume)
    search.update(source_id=state[0], metal=state[1], protocol_id=PROTOCOL,
                  previous_restricted_search=loaded['terminals'][state]['search'],
                  original_source=loaded['terminals'][state]['source']['prepared']['physical_atoms'],
                  qualification_result=cfg['qualification_result'], geometry_source_rebased=False)
    write(d/'SEARCH.json', search)
    row = {'source_id': state[0], 'metal': state[1], 'status': search['status'],
           'physical_state_admitted': search['physical_state_admitted'], 'terminal_point': search['terminal_point'],
           'search': pin(d/'SEARCH.json'), 'points': evaluator.engine.records, 'requests': evaluator.engine.requests,
           'actual_cost': evaluator.engine.cost, 'scout_result': cfg.get('scout_result') if stage == 'run' else None}
    write(d/'RESULT.json', row); row['result'] = pin(d/'RESULT.json')
    result.update(status=search['status'], rows=[row], actual_costs=[evaluator.engine.cost],
                  continuation_execution_and_physical_admission=search['status'] == 'complete' and search['physical_state_admitted'],
                  qualification_result=cfg['qualification_result'])
    if stage == 'run':
        result['scout_result'] = cfg['scout_result']
    require(search['distinct_evaluations'] <= LIMITS[stage]['evaluations']
            and search['accepted_steps'] <= LIMITS[stage]['accepted_steps'], 'Search budget exceeded')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', required=True, choices=('preflight', 'qualify', 'scout', 'run'))
    parser.add_argument('--state', choices=['_'.join(s) for s in STATES])
    parser.add_argument('--config', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) >= 1, 'Allocated CPU required')
    if args.stage != 'preflight':
        require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'), 'Allocated GPU required')
    require((args.state is not None) == (args.stage == 'run'), '--state is required only for run')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    cfg = json.loads(args.config.read_text()); tick = time.monotonic()
    result = {'protocol_id': PROTOCOL, 'stage': args.stage, 'status': 'execution_failure',
              'config': pin(args.config), 'implementation': pin(__file__), 'plan': cfg['plan'], 'spec': cfg['spec'],
              'optimizer': cfg['optimizer'], 'parent_pool_result': cfg['parent_pool_result'],
              'admission_revision':'archived_metric_parent_replay_v2',
              'replay_helper':cfg['replay_helper'],'recovery_evidence':cfg['recovery_evidence'],
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'rows': [], 'actual_costs': [],
              'new_GFN_calls': 0, 'new_MM_calls': 0, 'new_DFT_calls': 0, 'new_initialization_calls': 0}
    write(out/'START.json', result)
    try:
        loaded = load_inputs(cfg, out)
        result.update(preflight=loaded['preflight'], settings=loaded['core'].SETTINGS,
                      geometry_helper=loaded['geometry_helper'],
                      archived_metric_replay=loaded['archived_metric_replay'])
        write(out/'PREFLIGHT.json', loaded['preflight']); result['preflight_result'] = pin(out/'PREFLIGHT.json')
        if args.stage == 'preflight':
            result['status'] = 'complete'
        elif args.stage == 'qualify':
            qualify(cfg, loaded, out, result)
        else:
            state = STATES[0] if args.stage == 'scout' else tuple(args.state.split('_'))
            run_search(cfg, loaded, out, result, state, args.stage)
        verify(result['config']); verify(cfg['plan']); verify(cfg['spec']); verify(cfg['optimizer']); verify(cfg['replay_helper']); verify(cfg['recovery_evidence'])
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick,
                  allocation={k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'CUDA_VISIBLE_DEVICES')})
    result['new_composite_evaluations'] = sum(c['new_composite_evaluations'] for c in result['actual_costs'])
    write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
