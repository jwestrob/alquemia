#!/usr/bin/env python3
"""Source-specific full-Cartesian LBT common pools; at most four cross energies."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
import traceback

PROTOCOL = 'lbt_cartesian_source_pools_v1'
STATES = (('7CCO', 'La'), ('7CCO', 'Dy'), ('7CCN', 'La'), ('7CCN', 'Dy'))
METALS = ('La', 'Dy')
SOURCES = ('7CCO', '7CCN')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: '+record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def coordinate_key(sid, x):
    return sid, hashlib.sha256(np.asarray(x, dtype='<f8').tobytes(order='C')).hexdigest()


def menu(helper, sid, ids, cells, slots):
    missing = [{'candidate_id': cid, 'metal': metal} for cid in ids for metal in METALS
               if cells.get((cid, metal), {}).get('status') != 'complete']
    value = {'source_scope': sid, 'status': 'unavailable', 'candidate_ids': ids,
             'declared_states': 2*len(ids), 'declared_search_slots': len(slots), 'missing_cells': missing}
    if missing or any(s['status'] != 'admitted' for s in slots):
        value['reason'] = 'A required terminal or paired state is unavailable'
        return value
    selected = {m: min((cells[(cid, m)] for cid in ids), key=lambda r: r['total_energy_kcal_mol']) for m in METALS}
    native = {m: selected[m]['native_energies'] for m in METALS}
    diff = helper.difference(native['Dy'], native['La'])
    value.update(status='available', selected={m: selected[m]['candidate_id'] for m in METALS},
        exact_ties={m: [cid for cid in ids if cells[(cid, m)]['total_energy_kcal_mol']
                        == selected[m]['total_energy_kcal_mol']] for m in METALS},
        selected_native_energies=native, LBT_Dy_minus_La_components_kcal_mol=diff,
        LBT_Dy_minus_La_kcal_mol=math.fsum(diff.values()))
    return value


def replay_baseline(loaded, out):
    """All44 actual cells, including original energy-only cells; no new engines."""
    base, world = loaded['base'], loaded['world']; parent = loaded['parent_pool']; helper = world['h']
    candidates = copy.deepcopy(parent['candidates'])
    cells = {(r['candidate_id'], m): copy.deepcopy(c) for r in parent['paired_matrix'] for m, c in r['cells'].items()}
    require(len(candidates) == 22 and len(cells) == 44, 'Old22/44 denominator differs')
    by_coordinate = {}; receipts = []
    for candidate in candidates:
        sid = candidate['source_id']; source = world['sources'][sid]
        physical = read(candidate['physical_source']); x = np.asarray([a['xyz_A'] for a in physical])
        require(len(physical) == 256 and
                [{k: v for k, v in a.items() if k != 'xyz_A'} for a in physical] ==
                [{k: v for k, v in a.items() if k != 'xyz_A'} for a in source['atoms']],
                'Old candidate physical membership differs')
        key = coordinate_key(sid, x); require(key not in by_coordinate, 'Old exact-coordinate duplicate')
        candidate['_x'] = x; by_coordinate[key] = candidate
        for metal in METALS:
            row = cells[(candidate['candidate_id'], metal)]
            require(row['status'] == 'complete' and row['source_id'] == sid and row['metal'] == metal,
                    'Old paired row identity differs')
            if 'point' in row:
                _, rx, native = base.admit_point(world, source, metal, row['point'], False)
                require(np.array_equal(x, rx), 'Old POINT coordinates differ')
                evidence = row['point']
            else:
                ev = row['MACE_evidence']; raw = read(ev['raw_MACE'])
                c = world['matched_cells'][(candidate['candidate_id'], metal)]
                if ev['kind'] == 'origin':
                    raw = next(r for r in raw['rows'] if r['task_id'] == c['task']['task_id'])
                d = out/(candidate['candidate_id']+'_'+metal); d.mkdir()
                write(d/'MACE.json', raw)
                point = {'MACE': pin(d/'MACE.json'), 'ALPB': row['ALPB'],
                         'xyz': ev['xyz'], 'physical_source': ev['physical_source']}
                require(ev['physical_source'] == candidate['physical_source'] and np.array_equal(x, c['x']),
                        'Old evidence physical source differs')
                native, _ = base.raw_components(world, source, metal, point, False)
                evidence = ev['evidence_record']
            parts = helper.native_components(native)
            require(native == row['native_energies'] and parts == row['component_energies_kcal_mol']
                    and math.fsum(parts.values()) == row['total_energy_kcal_mol'], 'Old native energy replay differs')
            receipts.append({'candidate_id': candidate['candidate_id'], 'metal': metal, 'evidence': evidence, 'pass': True})
    baselines = []
    for sid in SOURCES:
        saved = next(r for r in parent['comparison']['expanded_menus'] if r['source_scope'] == sid)
        ids = [c['candidate_id'] for c in candidates if c['source_id'] == sid]
        require(ids == saved['candidate_ids'], 'Old source menu order differs')
        replay = menu(helper, sid, ids, cells, [])
        require(replay['status'] == 'available' and all(replay[k] == saved[k] for k in
                ('selected', 'exact_ties', 'selected_native_energies', 'LBT_Dy_minus_La_components_kcal_mol',
                 'LBT_Dy_minus_La_kcal_mol')), 'Old source menu selection/contrast differs')
        baselines.append(copy.deepcopy(saved))
    # Real saved-artifact test of the same admission routine used for new menus.
    sid = SOURCES[0]; ids = baselines[0]['candidate_ids']; removed = (ids[0], 'La')
    partial = dict(cells); partial.pop(removed)
    missing_test = menu(helper, sid, ids, partial, [])
    require(missing_test['status'] == 'unavailable' and
            {'candidate_id': removed[0], 'metal': removed[1]} in missing_test['missing_cells'],
            'Real missing-required-cell test failed')
    check = {'status': 'passed', 'new_molecular_calls': 0, 'old_geometries': 22, 'old_cells': 44,
             'cell_replays': receipts, 'old_source_menus': baselines,
             'missing_cell_test': {'test_only': True, 'removed_actual_cell': list(removed), 'result': missing_test}}
    write(out/'BASELINE_REPLAY.json', check)
    return candidates, cells, by_coordinate, baselines, pin(out/'BASELINE_REPLAY.json')


def admit_runs(cfg, cart, ccfg, loaded, out, candidates, cells, by_coordinate):
    require(len(cfg['run_results']) == 4, 'Four ordered run result pins required')
    base, world = loaded['base'], loaded['world']; helper = world['h']; slots = []; records = []; requests = []
    for state, parent_pin in zip(STATES, cfg['run_results']):
        sid, metal = state; parent = read(parent_pin); pcfg = read(parent['config'])
        verify(parent['implementation']); verify(pcfg['implementation'])
        require(parent['implementation']['sha256'] == pcfg['implementation']['sha256'] == cfg['cartesian_implementation']['sha256']
                and parent['protocol_id'] == cart.PROTOCOL and parent['stage'] == 'run', 'Cartesian run lineage differs')
        keys = ['plan', 'spec', 'optimizer', 'parent_pool_result']
        keys += [k for k in ('replay_helper', 'recovery_evidence') if k in ccfg]
        for k in keys:
            verify(pcfg[k]); require(pcfg[k]['sha256'] == ccfg[k]['sha256'], 'Run recipe differs: '+k)
        slot = {'source_id': sid, 'metal': metal, 'run_result': parent_pin, 'status': 'unavailable'}; slots.append(slot)
        source = world['sources'][sid]
        records.extend((state, p) for p in loaded['terminals'][state]['records'])
        requests.extend((state, r) for r in loaded['terminals'][state]['owner']['requests'])
        qualification = cart.parent_gate(pcfg, 'qualification_result', 'qualify')
        require(qualification['qualification_passed'], 'Run qualification unavailable')
        scout = cart.parent_gate(pcfg, 'scout_result', 'scout')
        require(scout['continuation_execution_and_physical_admission'] and len(scout['rows']) == 1,
                'Original scout admission differs')
        if state == STATES[0]:
            records.extend((state, p) for p in scout['rows'][0]['points'])
            requests.extend((state, r) for r in scout['rows'][0]['requests'])
        if len(parent['rows']) != 1:
            slot['reason'] = 'Run did not publish exactly one terminal row'; continue
        row = parent['rows'][0]
        require((row['source_id'], row['metal']) == state and read(row['result']) ==
                {k: v for k, v in row.items() if k != 'result'}, 'Run row identity/receipt differs')
        records.extend((state, p) for p in row['points']); requests.extend((state, r) for r in row['requests'])
        search = read(row['search']); slot.update(search=row['search'], termination=search['termination'],
            full_force=search['final_full_force'], accepted_steps=search['accepted_steps'],
            distinct_evaluations=search['distinct_evaluations'], converged=search['converged'],
            work_from_start=search['work_from_start'], final_geometry=search['final_geometry'])
        checkpoint = read(search['checkpoint'])
        require(checkpoint['point'] == search['terminal_point'] and checkpoint['initial_point'] == search['start_point']
                and checkpoint['accepted_steps'] == search['accepted_steps']
                and checkpoint['distinct_evaluations'] == search['distinct_evaluations']
                and checkpoint['settings'] == search['settings'] and checkpoint['limits'] == search['limits'],
                'Final cumulative checkpoint differs')
        records.extend((state, p) for p in checkpoint['evaluated_points'])
        requests.extend((state, r) for r in checkpoint['evaluation_requests'])
        require(search['start_point'] == loaded['terminals'][state]['point_pin']
                and search['original_source'] == source['prepared']['physical_atoms']
                and search['previous_restricted_search'] == loaded['terminals'][state]['search']
                and search['geometry_source_rebased'] is False and search['limits'] == cart.LIMITS['run']
                and search['settings'] == loaded['core'].SETTINGS and search['accepted_steps'] <= 80
                and search['distinct_evaluations'] <= 200, 'Search original source/budget/settings differ')
        require((search['source_id'], search['metal']) == state and search['protocol_id'] == cart.PROTOCOL
                and row['terminal_point'] == search['terminal_point'], 'Search endpoint identity differs')
        if parent['status'] != 'complete' or row['status'] != 'complete' or search['status'] != 'complete' or not (
                row['physical_state_admitted'] and search['physical_state_admitted']):
            slot['reason'] = 'Execution or original bound-state admission failed'; continue
        point, x, native = base.admit_point(world, source, metal, search['terminal_point'], True)
        check = loaded['guards'][state].check(x)
        gradient = np.load(verify(point['gradient_kcal_mol_A']), allow_pickle=False)
        require(check == search['final_geometry'] and check['pass'] and check['state_compatible']
                and loaded['core'].force_metrics(gradient) == search['final_full_force'], 'Terminal guard/full-force replay differs')
        key = coordinate_key(sid, x); candidate = by_coordinate.get(key)
        if candidate is None:
            candidate = {'candidate_id': 'cartesian_%02d' % len(slots), 'source_id': sid,
                         'physical_source': point['physical_source'], 'physical_coordinate_bytes_sha256': key[1],
                         'aliases': [], '_x': x}
            candidates.append(candidate); by_coordinate[key] = candidate
        require(coordinate_key(sid, candidate['_x']) == key and np.array_equal(candidate['_x'], x), 'Coordinate hash collision')
        alias = sid+':cartesian_terminal_'+metal
        require(alias not in candidate['aliases'], 'Duplicate new alias'); candidate['aliases'].append(alias)
        cell_key = (candidate['candidate_id'], metal)
        if cell_key in cells:
            require(cells[cell_key]['native_energies'] == native, 'Conflicting exact own-state components')
        else:
            parts = helper.native_components(native)
            cells[cell_key] = {'candidate_id': candidate['candidate_id'], 'source_id': sid, 'metal': metal,
                'status': 'complete', 'point': search['terminal_point'], 'native_energies': native,
                'component_energies_kcal_mol': parts, 'total_energy_kcal_mol': math.fsum(parts.values()), 'reused': True}
        slot.update(status='admitted', point=search['terminal_point'], candidate_id=candidate['candidate_id'], alias=alias)
    require(len(candidates) <= 26, 'Four-new-geometry limit exceeded')
    seen = set()
    for state, record in records:
        identity = (state, record['sha256'])
        if identity in seen:
            continue
        seen.add(identity); sid, metal = state
        p, x, native = base.admit_point(world, world['sources'][sid], metal, record, False)
        candidate = by_coordinate.get(coordinate_key(sid, x))
        if candidate is None:
            continue  # A trial can supply a cell, never an extra candidate geometry.
        cid = candidate['candidate_id']; key = (cid, metal)
        if key in cells:
            require(cells[key]['native_energies'] == native, 'Conflicting exact cached component')
        else:
            parts = helper.native_components(native)
            cells[key] = {'candidate_id': cid, 'source_id': sid, 'metal': metal, 'status': 'complete',
                'point': record, 'native_energies': native, 'component_energies_kcal_mol': parts,
                'total_energy_kcal_mol': math.fsum(parts.values()), 'reused': True, 'trial_not_added_as_candidate': True}
    blocked = {}; seen = set()
    for state, request in requests:
        rp = request['record']
        if rp['sha256'] in seen:
            continue
        seen.add(rp['sha256']); sid, metal = state; source = world['sources'][sid]
        require(read(rp) == {k: v for k, v in request.items() if k != 'record'} and
                request['source_id'] == sid and request['metal'] == metal and
                request['recipe'] == base.recipe(world, source, metal), 'Saved request identity differs')
        x = world['pool'].check_coordinates(world['kinworld'], source, metal,
                                            request['xyz'], request['physical_source'], request['positions'])
        require(request['cache_key'] == base.key_for(world, source, metal, x), 'Saved request coordinate key differs')
        if request['status'] == 'complete':
            base.admit_point(world, source, metal, request['point'], False); continue
        require(request['status'] == 'unavailable', 'Unknown request outcome')
        candidate = by_coordinate.get(coordinate_key(sid, x))
        if candidate is not None and (candidate['candidate_id'], metal) not in cells:
            blocked[(candidate['candidate_id'], metal)] = {'request': rp, 'reason': 'Earlier exact-state failure; no retry'}
    return slots, blocked


def comparison(helper, candidates, cells, slots, baselines):
    menus = []
    for sid in SOURCES:
        required = [s for s in slots if s['source_id'] == sid]
        ids = [c['candidate_id'] for c in candidates if c['source_id'] == sid]
        value = menu(helper, sid, ids, cells, required); old = next(b for b in baselines if b['source_scope'] == sid)
        value['old_baseline'] = old
        if value['status'] == 'available':
            value['selected_work_from_old_minimum_components_kcal_mol'] = {
                m: helper.difference(value['selected_native_energies'][m], old['selected_native_energies'][m]) for m in METALS}
            value['selected_work_from_old_minimum_kcal_mol'] = {
                m: math.fsum(v.values()) for m, v in value['selected_work_from_old_minimum_components_kcal_mol'].items()}
            value['contrast_change_from_old_kcal_mol'] = math.fsum(
                value['LBT_Dy_minus_La_components_kcal_mol'][k]-old['LBT_Dy_minus_La_components_kcal_mol'][k]
                for k in helper.NAMES)
        menus.append(value)
    old_signed = math.fsum((baselines[1]['LBT_Dy_minus_La_kcal_mol'], -baselines[0]['LBT_Dy_minus_La_kcal_mol']))
    spread = {'status': 'unavailable', 'old_7CCN_minus_7CCO_kcal_mol': old_signed}
    if all(m['status'] == 'available' for m in menus):
        signed = math.fsum((menus[1]['LBT_Dy_minus_La_kcal_mol'], -menus[0]['LBT_Dy_minus_La_kcal_mol']))
        spread.update(status='available', new_7CCN_minus_7CCO_kcal_mol=signed, absolute_spread_kcal_mol=abs(signed))
    return {'source_menus': menus, 'source_spread': spread, 'no_cross_source_selection': True,
            'interpretation': 'Conditional finite bound-state energy menus; residual forces retained; no affinity or populations.'}


def run(cfg, out, result, stage):
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Pool implementation differs')
    verify(cfg['plan']); verify(cfg['policy']); cart = module('cartesian_pool_adapter', cfg['cartesian_implementation'])
    ccfg = read(cfg['cartesian_config']); verify(ccfg['implementation'])
    require(ccfg['implementation']['sha256'] == cfg['cartesian_implementation']['sha256'] and
            ccfg['plan']['sha256'] == cfg['plan']['sha256'], 'Cartesian source/plan differs')
    d = out/'source_admission'; d.mkdir(); loaded = cart.load_inputs(ccfg, d)
    global np
    import numpy as np
    d = out/'old_baseline'; d.mkdir()
    candidates, cells, by_coordinate, baselines, replay = replay_baseline(loaded, d)
    old_cells = copy.deepcopy(cells); result.update(baseline_replay=replay, old_source_baselines=baselines,
        old_parent_pool=ccfg['parent_pool_result'], old_geometries=22, old_cells=44,
        actual_cost={k: 0 for k in loaded['base'].COST_KEYS})
    if stage == 'baseline':
        result.update(status='complete', declared_geometries=22, declared_cells=44,
                      model_initialized=False, baseline_only=True)
        return
    slots, blocked = admit_runs(cfg, cart, ccfg, loaded, out, candidates, cells, by_coordinate)
    helper, base, world = loaded['world']['h'], loaded['base'], loaded['world']
    missing = [(c, m) for c in candidates for m in METALS if (c['candidate_id'], m) not in cells
               and (c['candidate_id'], m) not in blocked]
    require(len(missing) <= 4 and all(c['candidate_id'].startswith('cartesian_') for c, _ in missing),
            'Only four missing new-terminal opposite-state cells are authorized')
    result.update(terminal_slots=slots, blocked_exact_failures=[dict(candidate_id=cid, metal=m, **v) for (cid, m), v in blocked.items()],
                  declared_geometries=len(candidates), declared_cells=2*len(candidates),
                  actual_cost={k: 0 for k in base.COST_KEYS})
    all_admitted = all(s['status'] == 'admitted' for s in slots)
    write(out/'PREFLIGHT.json', {'status': 'passed', 'old_baseline_replay': replay, 'new_molecular_calls': 0,
        'model_initialized': False, 'terminal_slots': slots, 'all_four_terminals_admitted': all_admitted,
        'required_cross_cells': [{'candidate_id': c['candidate_id'], 'metal': m} for c, m in missing],
        'blocked_exact_failures': result['blocked_exact_failures'], 'new_cross_ceiling': 4})
    result['preflight'] = pin(out/'PREFLIGHT.json'); native = None
    for c, metal in missing:
        row = {'candidate_id': c['candidate_id'], 'source_id': c['source_id'], 'metal': metal,
               'status': 'unavailable', 'reused': False}; result['cross_rows'].append(row)
        if not all_admitted:
            row['reason'] = 'No new calls before all four required terminal admissions'; continue
        require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'), 'Allocated GPU required for cross cells')
        d = out/(c['candidate_id']+'_'+metal); d.mkdir(); source = world['sources'][c['source_id']]
        evaluator = base.MatchedEvaluator(world, source, metal, d, native=native)
        try:
            check = loaded['guards'][(c['source_id'], metal)].check(c['_x'])
            require(check['pass'] and check['state_compatible'], 'Cross candidate physical/state gate failed')
            result['new_cross_cells_attempted'] += 1
            point, record, hit = evaluator.evaluate(c['_x'], compute_force=False)
            require(not hit, 'Undeclared hidden cross cache')
            point, x, values = base.admit_point(world, source, metal, record, False)
            require(coordinate_key(c['source_id'], x) == coordinate_key(c['source_id'], c['_x']), 'Cross coordinates differ')
            raw = read(point['MACE']); require(raw['forces_eV_A'] is None and raw['force_convention'] is None, 'Unexpected cross force request')
            parts = helper.native_components(values)
            row.update(status='complete', point=record, native_energies=values,
                       component_energies_kcal_mol=parts, total_energy_kcal_mol=math.fsum(parts.values()))
            cells[(c['candidate_id'], metal)] = row
        except Exception as exc:
            row.update(failure=repr(exc), traceback=traceback.format_exc())
        finally:
            native = evaluator.native; row.update(actual_cost=dict(evaluator.cost), requests=evaluator.requests)
            for k in base.COST_KEYS:
                result['actual_cost'][k] += evaluator.cost[k]
    require(result['new_cross_cells_attempted'] <= 4 and all(result['actual_cost'][k] <= 4 for k in
            ('new_composite_evaluations', 'MACE_calls_attempted', 'ALPB_calls_attempted')), 'Cross call ceiling exceeded')
    require(all(cells[k] == v for k, v in old_cells.items()), 'An old cell changed')
    result['comparison'] = comparison(helper, candidates, cells, slots, baselines)
    result['paired_matrix'] = [{'candidate_id': c['candidate_id'], 'source_id': c['source_id'], 'aliases': c['aliases'],
        'cells': {m: cells.get((c['candidate_id'], m), {'status': 'unavailable'}) for m in METALS}} for c in candidates]
    result['candidates'] = [{k: v for k, v in c.items() if k != '_x'} for c in candidates]
    result['status'] = 'complete' if all_admitted and all(m['status'] == 'available'
                        for m in result['comparison']['source_menus']) else 'unavailable'
    if native is not None:
        require(native.calls == result['actual_cost']['MACE_calls_attempted'], 'Warm native call count differs')
        result.update(native_model_calls=native.calls, native_device=native.device)
    lines = ['# LBT full-coordinate source pools', '', 'Status: '+result['status']+'.', '',
             '| Source | Availability | La work | Dy work | Contrast change |', '| --- | --- | ---: | ---: | ---: |']
    for m in result['comparison']['source_menus']:
        if m['status'] == 'available':
            w = m['selected_work_from_old_minimum_kcal_mol']
            lines.append('| %s | available | %.6f | %.6f | %.6f |' %
                         (m['source_scope'], w['La'], w['Dy'], m['contrast_change_from_old_kcal_mol']))
        else:
            lines.append('| '+m['source_scope']+' | unavailable | — | — | — |')
    lines += ['', 'Works and contrast changes are kcal/mol. Original source baselines, all failures, raw full forces,',
              'termination and exact ties are retained in RESULT.json. These are conditional finite energies, not affinities.', '']
    (out/'REPORT.md').write_text('\n'.join(lines)); result['report'] = pin(out/'REPORT.md')


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--stage', choices=('baseline', 'run'), default='run'); args = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); cfg = json.loads(args.config.read_text()); tick = time.monotonic()
    result = {'protocol_id': PROTOCOL, 'stage': args.stage, 'status': 'execution_failure', 'implementation': pin(__file__), 'config': pin(args.config),
              'plan': cfg['plan'], 'policy': cfg['policy'], 'run_results': cfg.get('run_results'),
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'declared_terminal_slots': 0 if args.stage == 'baseline' else 4,
              'new_cross_cells_attempted': 0, 'cross_rows': [], 'new_GFN_calls': 0, 'new_MM_calls': 0,
              'new_DFT_calls': 0, 'optimizer_steps': 0, 'new_MACE_gradients_requested': 0}
    try:
        run(cfg, out, result, args.stage); verify(result['config']); verify(cfg['policy'])
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic()-tick; write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
