#!/usr/bin/env python3
"""One explicitly chosen LBT hybrid region, unchanged Cartesian search and paired pool."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import time
import traceback

PROTOCOL = 'lbt_one_region_full_physical_hybrid_response_v1'
STATES = (('7CCO', 'La'), ('7CCO', 'Dy'), ('7CCN', 'La'), ('7CCN', 'Dy'))
LIMITS = {'scout': {'accepted_steps': 12, 'evaluations': 40},
          'run': {'accepted_steps': 80, 'evaluations': 200}}
FROZEN = {
    'optimizer': 'd93821da485f8d078e42bf7031018931632858ea720af617f26da2dab77617ee',
    'geometry_helper': '07f29ca60adefa185ad89ec68a51fe3221149230c80d016c00d6994ceec98226',
    'plan': '9e49c81ac4990c76f5bff6bd8c4e445029d3c93761761e9db4f14ec26b85fe03',
}
DONORS = ('A/3//OD1', 'A/3//OD2', 'A/5//OD1', 'A/7//OD1', 'A/7//OD2',
          'A/9//O', 'A/11//OE1', 'A/11//OE2', 'A/14//OE1', 'A/14//OE2')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        sha = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Changed pin: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
    return pin(path)


def module(name, p):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def same(a, b):
    verify(a)
    verify(b)
    return a['sha256'] == b['sha256']


def prepare(cfg, out, result):
    """Replay qualified chemistry and the complete 104-cell parent. Zero new calls."""
    require(cfg['protocol'] == PROTOCOL and same(cfg['implementation'], pin(__file__)), 'Response source/protocol differs')
    require(cfg['region'] in ('baseline', 'expanded') and cfg['states'] == [list(s) for s in STATES],
            'Exactly one region and all four states required')
    verify(cfg['specification'])
    for key, sha in FROZEN.items():
        require(verify(cfg[key]) and cfg[key]['sha256'] == sha, 'Frozen helper differs: ' + key)
    parent = read(cfg['parent_pool_result'])
    pcfg = read(parent['config'])
    require(parent['status'] == 'complete' and parent['stage'] == 'run'
            and parent['declared_geometries'] == 26 and parent['declared_hybrid_cells'] == 104
            and same(parent['implementation'], cfg['parent_pool_implementation'])
            and same(pcfg['implementation'], cfg['parent_pool_implementation'])
            and same(pcfg['plan'], cfg['plan']), 'Actual complete 26-geometry/104-cell hybrid pool required')
    pool = module('hybrid_response_parent_pool', cfg['parent_pool_implementation'])
    require(parent['protocol'] == pool.PROTOCOL, 'Parent pool protocol differs')
    scratch = {}
    directory = out / 'parent_replay'
    directory.mkdir()
    try:
        base, ev, candidates, cells, old_cells, whole = pool.prepare(pcfg, directory, scratch)
    finally:
        if '_evaluator' in scratch:
            result['_evaluator'] = scratch['_evaluator']
    qualifier = module('hybrid_response_restore', pcfg['qualification_implementation'])
    require(base.pin(base.__file__)['sha256'] == qualifier.EVALUATOR_SHA
            == '686eba84fd3dfcf28724d2d4e7957528c987445d7137b18c88a07dfec055f661', 'Frozen hybrid evaluator differs')
    require(parent['identity_sha256'] == ev.identity_hash and parent['candidates'] == scratch['candidates'],
            'Full parent candidate/model identity differs')
    expected = [(r['candidate_id'], r['region'], r['metal']) for r in scratch['cells']]
    actual = [(r['candidate_id'], r['region'], r['metal']) for r in parent['cells']]
    require(actual == expected and len(actual) == len(set(actual)) == 104, 'Complete ordered parent denominator differs')
    by_id = {c['candidate_id']: c for c in candidates}
    for slot, key in zip(parent['cells'], actual):
        require(slot['status'] == 'complete', 'Unavailable parent hybrid cell')
        p = qualifier.restore(ev, slot['point'])
        c = by_id[key[0]]
        require((p['source_id'], p['region'], p['metal']) == (c['source_id'], key[1], key[2])
                and ev.np.array_equal(p['positions_A'], c['_x'])
                and all(slot[k] == p[k] for k in ('native_energies', 'component_energies_kcal_mol', 'total_energy_kcal_mol'))
                and slot['cache_key'] == p['provenance']['cache_key'], 'Parent cell native state differs')
        cells[key] = copy.deepcopy(slot)
    require(pool.summarize(base, candidates, cells, old_cells, whole, scratch)
            and scratch['comparison'] == parent['comparison'], 'Parent selections/work do not replay')
    core = module('hybrid_response_cartesian_core', cfg['optimizer'])
    geometry = module('hybrid_response_geometry', cfg['geometry_helper'])
    region = cfg['region']
    starts, guards, metadata = {}, {}, []
    for sid, metal in STATES:
        menu = next(m for m in parent['comparison']['complete_source_menus']
                    if (m['source_id'], m['region']) == (sid, region))
        require(menu['status'] == 'available', 'Required source-local selected menu unavailable')
        cid = menu['selected'][metal]
        start = qualifier.restore(ev, cells[(cid, region, metal)]['point'])
        _, part, _ = ev.partitions[(sid, region)]
        original = read(part['source_preparation'])
        atoms = part['physical_atoms']
        ids = [a['id'] for a in atoms]
        require(len(atoms) == 256 and ids[-1] == 'metal' and original['protein_atoms'] == 255
                and original['waters'] == original['spectators'] == original['caps'] == [],
                'Whole LBT chemical inventory changed')
        source_atoms = read(original['physical_atoms'])
        require(source_atoms == atoms and ev.np.array_equal([a['xyz_A'] for a in source_atoms], ev.origins[(sid, region)]),
                'Original full physical source differs')
        bonds = read(original['covalent_bonds'])
        pairs = [b['native_particle_indices'] for b in bonds]
        for b, pair in zip(bonds, pairs):
            require(len(pair) == 2 and all(0 <= i < 255 for i in pair)
                    and [ids[i] for i in pair] == [b['atom_a_id'], b['atom_b_id']], 'Full covalent graph mapping differs')
        donor_map = read(original['endpoints']['matched_core'][metal]['source_cap_map'])
        donor_ids = [f"{a['chain']}/{a['resnum']}//{a['atom']}" for a in donor_map['contacts']]
        require(len(donor_ids) == 10 and set(donor_ids) == set(DONORS), 'Declared site donor inventory differs')
        target = {'index': 255, 'label': 'LBT_target', 'site_donor_indices': [ids.index(x) for x in DONORS],
                  'coordinating_water_indices': []}
        guard = core.PhysicalGuard(atoms, pairs, [b['forcefield_equilibrium_A'] for b in bonds],
                                   ev.origins[(sid, region)], geometry, list(range(255)), [target])
        origin_check = guard.check(ev.origins[(sid, region)])
        check = guard.check(start['positions_A'])
        require(origin_check['pass'] and origin_check['state_compatible'] and check['pass'] and check['state_compatible'],
                'Original or selected start is physically unavailable')
        starts[(sid, metal)] = start
        guards[(sid, metal)] = guard
        metadata.append({'source_id': sid, 'metal': metal, 'region': region, 'selected_candidate_id': cid,
                         'exact_ties': menu['exact_ties'][metal], 'point': start['provenance'], 'geometry': check,
                         'initial_full_force': core.force_metrics(ev.np.asarray(start['gradient_kcal_mol_A'])),
                         'original_preparation': part['source_preparation'], 'full_covalent_graph': original['covalent_bonds'],
                         'waters': [], 'spectators': [], 'source_conditioning': 'Lu' if sid == '7CCN' else 'La'})
    require(ev.cost['composites_attempted'] == 0 and ev.native is None, 'Preparation made a molecular call')
    result.update(identity_sha256=ev.identity_hash, starts=metadata, parent_comparison=parent['comparison'],
                  parent_actual_cost=parent['actual_cost'], inherited_geometries=26, inherited_same_region_cells=52,
                  settings=core.SETTINGS)
    result['preflight'] = write(out / 'PREFLIGHT.json', {'status': 'passed', 'new_molecular_calls': 0,
        'parent': cfg['parent_pool_result'], 'identity_sha256': ev.identity_hash, 'region': region,
        'replayed_parent_cells': 104, 'starts': metadata, 'settings': core.SETTINGS,
        'all_atoms_movable': 256, 'new_directions': 0, 'source_motion_envelope': None})
    return {'base': base, 'ev': ev, 'pool': pool, 'qualifier': qualifier, 'core': core, 'geometry': geometry,
            'starts': starts, 'guards': guards, 'parent': parent, 'candidates': candidates, 'cells': cells}


class Evaluator:
    """Only state binding and saved-point admission; all physics stays in V2."""
    def __init__(self, loaded, state, region):
        self.loaded, self.state, self.region = loaded, state, region
        self.ev = loaded['ev']
        self.pin = pin

    @property
    def new_evaluations(self):
        return self.ev.cost['composites_attempted']

    @property
    def evaluation_requests(self):
        return self.ev.requests

    def read_point(self, provenance):
        p = self.loaded['qualifier'].restore(self.ev, provenance['point'])
        require(p['provenance'] == provenance and (p['source_id'], p['metal']) == self.state
                and p['region'] == self.region, 'Search point state/provenance differs')
        return p

    def evaluate(self, x):
        return self.ev.evaluate(self.state[0], self.region, self.state[1], x)

    def difference(self, new, old):
        return self.loaded['base'].difference(new, old)


def binding(cfg, loaded, state):
    return {'protocol': PROTOCOL, 'region': cfg['region'], 'source_id': state[0], 'metal': state[1],
            'identity_sha256': loaded['ev'].identity_hash,
            'pins': {k: cfg[k]['sha256'] for k in ('implementation', 'specification', 'plan', 'optimizer',
                                                'geometry_helper', 'parent_pool_result', 'parent_pool_implementation')},
            'initial_point_sha256': loaded['starts'][state]['provenance']['point']['sha256'],
            'settings': loaded['core'].SETTINGS}


def saved_search(record, cfg, loaded, state, stages):
    """Admit complete search evidence and every cached rejected trial; no retry."""
    r = read(record)
    rcfg = read(r['config'])
    require(r['protocol'] == PROTOCOL and r['stage'] in stages and r['status'] == 'complete'
            and r['continuation_execution_and_physical_admission']
            and same(r['implementation'], cfg['implementation'])
            and all(same(rcfg[k], cfg[k]) for k in ('implementation', 'specification', 'plan', 'optimizer',
                                                  'geometry_helper', 'parent_pool_result', 'parent_pool_implementation'))
            and rcfg['region'] == cfg['region'] and rcfg['states'] == cfg['states']
            and read(r['binding']) == binding(cfg, loaded, state), 'Search lineage or physical admission unavailable')
    search = read(r['search'])
    require((search['source_id'], search['metal']) == state and search['region'] == cfg['region']
            and search['status'] == 'complete' and search['physical_state_admitted']
            and search['failure'] is None and search['termination'] != 'execution_failure', 'Search state/failure differs')
    checkpoint = read(search['checkpoint'])
    require(checkpoint['failure'] is None and checkpoint['point'] == search['terminal_point']
            and checkpoint['initial_point'] == search['start_point'] == loaded['starts'][state]['provenance']
            and checkpoint['accepted_steps'] == search['accepted_steps']
            and checkpoint['distinct_evaluations'] == search['distinct_evaluations']
            and checkpoint['settings'] == search['settings'] == loaded['core'].SETTINGS
            and checkpoint['limits'] == search['limits'] == LIMITS[r['stage']], 'Final checkpoint differs')
    require(all(type(checkpoint[k]) is int and 0 <= checkpoint[k] <= LIMITS[r['stage']][limit]
                for k, limit in (('accepted_steps', 'accepted_steps'), ('distinct_evaluations', 'evaluations'))),
            'Saved finite search counters exceed declared limits')
    evaluator = Evaluator(loaded, state, cfg['region'])
    requests = checkpoint['evaluation_requests']
    require(len({p['sha256'] for p in requests}) == len(requests)
            and len(requests) == checkpoint['distinct_evaluations'], 'Cumulative request accounting differs')
    requested = set()
    for request in requests:
        raw = read(request)
        require(raw['status'] == 'complete', 'Recorded failed request cannot be retried')
        point = loaded['qualifier'].restore(loaded['ev'], request)
        evaluator.read_point(point['provenance'])
        requested.add(request['sha256'])
    for provenance in checkpoint['evaluated_points']:
        # A trial may hit an exact parent/scout cache entry and cost no new call.
        require(provenance['point']['sha256'] in requested or provenance['cache_key'] in loaded['ev'].cache,
                'Evaluated trial absent from requests and admitted parent cache')
        evaluator.read_point(provenance)
    for request in r['requests']:
        require(request['sha256'] in requested and read(request)['status'] == 'complete', 'Untracked allocation request')
    require(len(r['requests']) == r['actual_cost']['composites_attempted'], 'Allocation call accounting differs')
    evaluator.read_point(checkpoint['initial_point'])
    terminal = evaluator.read_point(checkpoint['point'])
    gate = loaded['guards'][state].check(terminal['positions_A'])
    require(gate == search['final_geometry'] and gate['pass'] and gate['state_compatible']
            and loaded['core'].force_metrics(loaded['ev'].np.asarray(terminal['gradient_kcal_mol_A'])) == search['final_full_force']
            and evaluator.difference(terminal, loaded['starts'][state]) == search['work_from_start'],
            'Terminal physical/force/work replay differs')
    return r, search, checkpoint, terminal


def run_search(cfg, loaded, out, result, state, stage):
    resume = None
    if stage == 'run':
        scout, _, cp, _ = saved_search(cfg['scout_result'], cfg, loaded, STATES[0], ('scout',))
        result['scout_result'] = cfg['scout_result']
        result['scout_cost_separate'] = scout['actual_cost']
        if state == STATES[0]:
            resume = cp
    if cfg.get('resume_result') is not None:
        require(stage == 'run', 'Explicit resume is permitted only for the unchanged run budget')
        _, _, resume, _ = saved_search(cfg['resume_result'], cfg, loaded, state, ('run',))
        result['resume_result'] = cfg['resume_result']
    result['binding'] = write(out / 'RUN_BINDING.json', binding(cfg, loaded, state))
    evaluator = Evaluator(loaded, state, cfg['region'])
    start = loaded['starts'][state]
    search = loaded['core'].search(evaluator, loaded['guards'][state], start, out / 'search', LIMITS[stage],
                                  loaded['geometry'].lbfgs_direction, resume=resume)
    search.update(source_id=state[0], metal=state[1], region=cfg['region'], protocol=PROTOCOL,
                  binding=result['binding'], start_policy='same_region_source_local_full_parent_pool_minimum')
    result['search'] = write(out / 'SEARCH.json', search)
    result.update(status=search['status'], source_id=state[0], metal=state[1],
                  continuation_execution_and_physical_admission=search['status'] == 'complete' and search['physical_state_admitted'],
                  terminal_point=search['terminal_point'], termination=search['termination'],
                  converged=search['converged'], final_full_force=search['final_full_force'],
                  accepted_steps=search['accepted_steps'], distinct_evaluations=search['distinct_evaluations'],
                  work_from_start=search['work_from_start'])
    require(search['accepted_steps'] <= LIMITS[stage]['accepted_steps']
            and search['distinct_evaluations'] <= LIMITS[stage]['evaluations'], 'Declared search budget exceeded')


def paired_pool(cfg, loaded, out, result):
    base, ev, pool = loaded['base'], loaded['ev'], loaded['pool']
    region = cfg['region']
    slots = [{'source_id': s, 'metal': m, 'status': 'unavailable', 'result': p}
             for (s, m), p in zip(STATES, cfg['run_results'])]
    require(len(cfg['run_results']) == len(slots) == 4, 'Four ordered search result slots required')
    result['terminal_slots'] = slots
    candidates = [{k: copy.deepcopy(v) for k, v in c.items()} for c in loaded['candidates']]
    cells = {k: copy.deepcopy(v) for k, v in loaded['cells'].items() if k[1] == region}
    by_coordinate = {(c['source_id'], ev.np.asarray(c['_x'], dtype='<f8').tobytes()): c for c in candidates}
    for slot, state in zip(slots, STATES):
        try:
            _, search, _, point = saved_search(slot['result'], cfg, loaded, state, ('run',))
            x = ev.np.asarray(point['positions_A'], dtype='<f8')
            key = (state[0], x.tobytes())
            alias = state[0] + ':hybrid_' + region + '_terminal_' + state[1]
            candidate = by_coordinate.get(key)
            if candidate is None:
                raw = read(point['provenance']['point'])
                physical = read(raw['task'])['physical_source']
                require(ev.np.array_equal([a['xyz_A'] for a in read(physical)], x), 'Terminal physical source differs')
                candidate = {'candidate_id': 'hybrid_terminal_%02d' % (STATES.index(state) + 1), 'source_id': state[0],
                             'physical_source': physical, 'physical_coordinate_bytes_sha256': hashlib.sha256(x.tobytes()).hexdigest(),
                             'aliases': [], '_x': x}
                candidates.append(candidate)
                by_coordinate[key] = candidate
            candidate['aliases'].append(alias)
            cid = candidate['candidate_id']
            own = {'candidate_id': cid, 'source_id': state[0], 'region': region, 'metal': state[1]}
            pool.cell_from_point(own, point, True)
            if (cid, region, state[1]) in cells:
                require(cells[(cid, region, state[1])]['native_energies'] == own['native_energies'], 'Duplicate coordinate energy differs')
            cells[(cid, region, state[1])] = own
            slot.update(status='admitted', candidate_id=cid, point=point['provenance'], search=read(slot['result'])['search'],
                        converged=search['converged'], full_force=search['final_full_force'], termination=search['termination'],
                        accepted_steps=search['accepted_steps'], distinct_evaluations=search['distinct_evaluations'])
        except Exception as exc:
            slot.update(failure=repr(exc), traceback=traceback.format_exc())
    result['candidates'] = [{k: v for k, v in c.items() if k != '_x'} for c in candidates]
    matrix = []
    for c in candidates:
        for metal in ('La', 'Dy'):
            key = (c['candidate_id'], region, metal)
            row = cells.setdefault(key, {'candidate_id': c['candidate_id'], 'source_id': c['source_id'], 'region': region,
                                         'metal': metal, 'status': 'unavailable_not_attempted', 'physical_source': c['physical_source']})
            cache_key = ev.key(c['source_id'], region, metal, c['_x'])
            if row['status'] != 'complete' and cache_key in ev.cache:
                require(ev.cache[cache_key] is not None, 'Failed exact cross cell cannot be retried')
                pool.cell_from_point(row, ev.cache[cache_key], True)
            matrix.append(row)
    result.update(paired_cells=matrix, declared_geometries=len(candidates), declared_cells=len(matrix), maximum_new_cross_composites=4)
    require(26 <= len(candidates) <= 30 and sum(r['status'] != 'complete' for r in matrix) <= 4,
            'Paired terminal pool exceeds four opposite-metal cells')
    if not all(s['status'] == 'admitted' for s in slots):
        result.update(status='unavailable_terminal_evidence', new_pool_available=False)
        return
    by_id = {c['candidate_id']: c for c in candidates}
    for row in matrix:
        if row['status'] == 'complete':
            continue
        before = len(ev.requests)
        row['status'] = 'unavailable_attempted'
        try:
            point = ev.evaluate(row['source_id'], region, row['metal'], by_id[row['candidate_id']]['_x'])
            pool.cell_from_point(row, point, False)
        except Exception as exc:
            row.update(failure=repr(exc), traceback=traceback.format_exc())
            break
        finally:
            row['attempted_point'] = ev.requests[-1] if len(ev.requests) > before else None
    menus = []
    for sid in ('7CCO', '7CCN'):
        old = next(m for m in loaded['parent']['comparison']['complete_source_menus']
                   if (m['source_id'], m['region']) == (sid, region))
        new = pool.menu(base, sid, region, [c['candidate_id'] for c in candidates if c['source_id'] == sid],
                        cells, loaded['parent']['origin_candidate_ids'][sid])
        row = {'source_id': sid, 'region': region, 'old': old, 'new': new}
        if new['status'] == 'available':
            works = {m: base.difference(cells[(new['selected'][m], region, m)], cells[(old['selected'][m], region, m)])
                     for m in ('La', 'Dy')}
            row.update(additional_selected_work_from_parent_minimum_kcal_mol=works,
                       change_in_Dy_minus_La_kcal_mol={k: works['Dy'][k] - works['La'][k] for k in works['La']})
        menus.append(row)
    result.update(comparison={'source_menus': menus, 'no_cross_source_pool': True, 'no_cross_region_pool': True},
                  new_pool_available=all(m['new']['status'] == 'available' for m in menus))
    result['status'] = 'complete' if result['new_pool_available'] else 'unavailable_backend'
    require(ev.cost['composites_attempted'] <= 4 and ev.cost['MACE_calls_attempted'] <= 4
            and ev.cost['GFN_calls_attempted'] <= 8 and ev.cost['MM_configurations_attempted'] <= 4
            and ev.cost['MM_group_queries_attempted'] <= 16 and ev.cost['ALPB_calls_attempted'] <= 4,
            'Four-cell cross-call component ceiling exceeded')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage', choices=('preflight', 'scout', 'run', 'pool'), required=True)
    p.add_argument('--state', choices=['_'.join(s) for s in STATES])
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required, including saved preflight')
    require((a.state is not None) == (a.stage == 'run'), '--state is required only for one run')
    if a.stage != 'preflight':
        require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'), 'Allocated GPU required')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick = time.monotonic()
    before = resource.getrusage(resource.RUSAGE_SELF)
    result = {'protocol': PROTOCOL, 'stage': a.stage, 'status': 'execution_failure', 'implementation': pin(__file__),
              'config': pin(a.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
              'declared_search_states': [list(s) for s in STATES], 'declared_search_denominator': 4,
              'new_DFT_calls': 0, 'new_derivative_directions': 0, 'unbound_draft_is_not_execution_authorization': True}
    try:
        cfg = json.loads(a.config.read_text())
        result.update(region=cfg['region'], parent_pool_result=cfg['parent_pool_result'])
        loaded = prepare(cfg, out, result)
        if a.stage == 'preflight':
            result['status'] = 'complete'
        elif a.stage == 'pool':
            paired_pool(cfg, loaded, out, result)
        else:
            state = STATES[0] if a.stage == 'scout' else tuple(a.state.split('_'))
            run_search(cfg, loaded, out, result, state, a.stage)
        verify(result['config'])
        for key in ('implementation', 'specification', 'plan', 'optimizer', 'geometry_helper', 'parent_pool_implementation', 'parent_pool_result'):
            verify(cfg[key])
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    ev = result.pop('_evaluator', None)
    if ev is not None:
        result.update(actual_cost=ev.cost, requests=ev.requests, points=ev.points)
        if ev.native is not None:
            result.update(native_MACE_calls=ev.native.calls, native_MACE_device=ev.native.device, native_MACE_versions=ev.native.versions)
    after = resource.getrusage(resource.RUSAGE_SELF)
    result.update(elapsed_seconds=time.monotonic() - tick,
                  parent_CPU_seconds=(after.ru_utime + after.ru_stime) - (before.ru_utime + before.ru_stime),
                  allocation={k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'CUDA_VISIBLE_DEVICES')})
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
