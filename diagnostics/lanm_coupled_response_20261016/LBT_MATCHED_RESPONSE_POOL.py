#!/usr/bin/env python3
"""Extend the frozen matched-FF LBT pool with four declared response terminals.

Source-only adapter; root owns execution. Old18/36 evidence and all15 old menus
remain explicit. At most four new opposite-state energy-only composite calls.
"""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import socket
import struct
import time
import traceback

SOURCES = ('7CCO', '7CCN')
METALS = ('La', 'Dy')
PROTOCOL = 'lbt3_matched_fixedFF_response_pool_v1'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def read(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: ' + record['path'])
    return json.loads(Path(record['path']).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def coordinate_hash(positions):
    return hashlib.sha256(b''.join(struct.pack('<d', 0. if v == 0. else v)
                                  for row in positions for v in row)).hexdigest()


def inputs(cfg, out):
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Pool implementation differs')
    require(pin(cfg['plan']['path']) == cfg['plan'], 'Plan differs')
    helper = cfg['response_helper']
    require(pin(helper['path'])['sha256'] == helper['sha256'], 'Response helper differs')
    spec = importlib.util.spec_from_file_location('matched_response_pool_helper', helper['path'])
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    parents = [read(cfg[key]) for key in ('scout_result', 'remaining_result')]
    configs = [read(parent['config']) for parent in parents]
    for parent, pcfg, stage in zip(parents, configs, ('scout', 'remaining')):
        require(parent['status'] in ('complete', 'unavailable') and parent['stage'] == stage,
                'Execution failure cannot supply terminal candidates')
        require(parent['implementation']['sha256'] == pcfg['implementation']['sha256'] == helper['sha256']
                and parent['plan']['sha256'] == pcfg['plan']['sha256'] == cfg['plan']['sha256']
                and parent['matched_pool_result'] == pcfg['matched_pool_result'] == cfg['matched_pool_result'],
                'Parent implementation/plan/model lineage differs')
        for record in (parent['implementation'], pcfg['implementation'], parent['plan'], pcfg['plan']):
            base.verify(record)
        pre = read(parent['preflight_result'])
        require(parent['protocol_id'] == base.PROTOCOL and pre == parent['preflight']
                and pre['status'] == 'passed' and pre['declared_states'] == 4
                and pre['molecular_calls'] == 0 and pre['model_loaded'] is False
                and len(pre['rows']) == 4 and all(r['pass'] for r in pre['rows'])
                and {(r['source_id'], r['metal']) for r in pre['rows']}
                == {(s, m) for s in SOURCES for m in METALS}, 'Parent four-start preflight differs')
        require((parent['status'] == 'complete')
                == all(r['status'] == 'complete' and r['candidate_admitted'] for r in parent['rows']),
                'Parent/row availability differs')
    require(parents[0]['settings'] == parents[1]['settings'], 'Parent search settings differ')
    require(parents[1]['scout_result'] == configs[1]['scout_result'] == cfg['scout_result'], 'Remaining scout gate differs')
    slots = [dict(row, parent=cfg[key]) for key, parent in zip(('scout_result', 'remaining_result'), parents)
             for row in parent['rows']]
    require(len(parents[0]['rows']) == 1 and (slots[0]['source_id'], slots[0]['metal']) == ('7CCO', 'La')
            and len(slots) == 4 and {(r['source_id'], r['metal']) for r in slots}
            == {(s, m) for s in SOURCES for m in METALS}, 'Four-search denominator differs')
    directory = out / 'source_admission'
    directory.mkdir()
    world = base.load_inputs(configs[0], directory)
    require(parents[0]['settings'] == world['settings']
            and parents[0]['fallback_rule'] == parents[1]['fallback_rule'] == base.FALLBACK,
            'Frozen search settings/fallback differ')
    prior = read(cfg['matched_pool_result'])
    require(world['matched_pool'] == prior and prior['status'] == 'complete'
            and prior['slurm_job_id'] == '12720' and prior['declared_geometries'] == 18
            and prior['declared_cells'] == len(prior['rows']) == 36
            and all(r['status'] == 'complete' for r in prior['rows']), 'Frozen18/36 pool differs')
    inventory = read(read(prior['config'])['inventory'])
    candidates = copy.deepcopy(inventory['candidates'])
    require(len(candidates) == 18 and [c['candidate_id'] for c in candidates]
            == ['geometry_%02d' % i for i in range(18)], 'Old candidate order differs')
    cells = {(r['candidate_id'], r['metal']): copy.deepcopy(r) for r in prior['rows']}
    require(len(cells) == 36 and set(cells) == {(c['candidate_id'], m) for c in candidates for m in METALS},
            'Old paired matrix differs')
    return base, world, prior, candidates, cells, slots, parents


def score_row(helper, candidate_id, source_id, metal, record, point, native):
    components = helper.native_components(native)
    return {'candidate_id': candidate_id, 'source_id': source_id, 'metal': metal, 'status': 'complete',
            'point': record, 'ALPB': point['ALPB'], 'native_energies': native,
            'component_energies_kcal_mol': components, 'total_energy_kcal_mol': math.fsum(components.values())}


def assemble(base, world, candidates, cells, slots):
    import numpy as np
    helper = world['h']
    by_coordinate = {}
    for candidate in candidates:
        physical = read(candidate['physical_source'])
        x = helper.vectors([a['xyz_A'] for a in physical], 256)
        require(coordinate_hash(x) == candidate['coordinate_sha256'], 'Old physical coordinate identity differs')
        candidate['_x'] = x
        by_coordinate[candidate['coordinate_sha256']] = candidate
    require(len(by_coordinate) == 18, 'Old exact-coordinate denominator differs')
    admitted = []
    start_ids = {('7CCO', 'La'): 'geometry_11', ('7CCO', 'Dy'): 'geometry_13',
                 ('7CCN', 'La'): 'geometry_15', ('7CCN', 'Dy'): 'geometry_17'}
    for row in sorted(slots, key=lambda r: (SOURCES.index(r['source_id']), METALS.index(r['metal']))):
        sid, metal = row['source_id'], row['metal']
        source = world['sources'][sid]
        start = world['starts'][(sid, metal)]
        require(read(row['result']) == {k: v for k, v in row.items() if k not in ('result', 'parent')},
                'Search slot receipt differs')
        require(all(row['actual_cost'][k] <= 24 for k in
                    ('new_composite_evaluations', 'MACE_calls_attempted', 'ALPB_calls_attempted')),
                'Per-search call budget differs')
        item = {'source_id': sid, 'metal': metal, 'parent': row['parent'], 'status': 'unavailable'}
        admitted.append(item)
        if row.get('search') is None:
            item['reason'] = 'Declared search has no terminal receipt'
            continue
        search = read(row['search'])
        item['search'] = row['search']
        require((search['source_id'], search['metal']) == (sid, metal), 'Search identity differs')
        if not row['candidate_admitted'] or not search['candidate_admitted']:
            item['reason'] = 'Declared search unavailable; trial points do not replace its terminal'
            continue
        require(row['status'] == 'complete' and search['protocol_id'] == base.PROTOCOL
                and search['status'] == 'finite_proposal_available'
                and search['candidate'] in row['points']
                and search['candidate'] == search['trace'][-1]['point']
                and len(search['trace']) == search['accepted_steps'] + 1
                and search['accepted_steps'] <= 12 and search['final_geometry']['pass']
                and search['initial_history_length'] == search['initial_new_evaluations'] == 0
                and search['initial_offset_v_A'] == [0.] * 10 and search['source_rebased'] is False
                and search['previous_search'] == row['previous_search'] == start['search_pin']
                and search['continuation_start'] == row['continuation_start'] == search['trace'][0]['point']
                and search['old_objective_terminal'] == row['old_objective_terminal'] == start['old_objective_terminal']
                and search['authoritative_origin'] == start['search']['authoritative_origin']
                and search['origin_guard']['pass'] and search['continuation_start_guard']['pass']
                and search['fallback_rule'] == base.FALLBACK
                and all(t['geometry']['pass'] and t['point'] in row['points'] for t in search['trace']),
                'Terminal search admission differs')
        first, first_x, first_native = base.admit_point(world, source, metal, search['trace'][0]['point'])
        old_id = start_ids[(sid, metal)]
        old_candidate = next(c for c in candidates if c['candidate_id'] == old_id)
        require(np.array_equal(first_x, np.asarray(old_candidate['_x']))
                and first_native == cells[(old_id, metal)]['native_energies'], 'Search matched start differs')
        point, x, native = base.admit_point(world, source, metal, search['candidate'])
        work = base.relative(point['components'], first['components'])
        require(work == search['work_from_continuation_start'] and work['composite_kcal_mol'] <= 1e-8,
                'Matched-start work replay differs or increased')
        q = np.asarray(search['full_q'], dtype=float)
        v = np.asarray(search['offset_v_A'], dtype=float)
        active = start['active']
        transform = start['transform']
        require(q.shape == (64,) and v.shape == (10,) and len(active) == 10
                and np.isfinite(q).all() and np.isfinite(v).all()
                and search['active_indices'] == list(active)
                and np.array_equal(np.asarray(search['coordinate_scaling_T']), transform)
                and np.array_equal(np.asarray(search['coordinate_offset_q64']), start['q'])
                and np.array_equal(q[[i for i in range(64) if i not in active]],
                                   start['q'][[i for i in range(64) if i not in active]])
                and np.allclose(q[active], start['q'][active] + transform @ v, rtol=0, atol=1e-12),
                'Original-source coordinate map/metric differs')
        replay, _, _, _ = source['kin'].evaluate(q)
        check = start['guard'].check(q, x, x)
        require(np.allclose(replay, x, rtol=0, atol=1e-12) and check['pass'], 'Terminal kinematic/physical replay failed')
        key = coordinate_hash(x)
        if key not in by_coordinate:
            candidate = {'candidate_id': 'geometry_%02d' % len(candidates), 'coordinate_sha256': key,
                         'source_id': sid, 'physical_source': point['physical_source'], 'aliases': [], '_x': x.tolist()}
            candidates.append(candidate)
            by_coordinate[key] = candidate
        candidate = by_coordinate[key]
        require(np.array_equal(np.asarray(candidate['_x']), x), 'Coordinate hash collision')
        alias = sid + ':matched_response_terminal_' + metal
        require(alias not in candidate['aliases'], 'Duplicate terminal alias')
        candidate['aliases'].append(alias)
        cell_key = (candidate['candidate_id'], metal)
        if cell_key in cells:
            require(cells[cell_key]['native_energies'] == native, 'Conflicting exact-state cached components')
        else:
            cells[cell_key] = score_row(helper, candidate['candidate_id'], sid, metal, search['candidate'], point, native)
            cells[cell_key]['reused'] = True
        item.update(status='admitted', alias=alias, candidate_id=candidate['candidate_id'],
                    point=search['candidate'], physical_check=check, accepted_steps=search['accepted_steps'],
                    termination=search['termination'], final_force_diagnostic=search['final_force_diagnostic'],
                    own_work_from_matched_start_components_kcal_mol=helper.difference(native, first_native))
    # A complete saved trial may supply a score at an already admitted candidate.
    # It never adds a geometry to the menu.
    seen = set()
    for row in slots:
        source = world['sources'][row['source_id']]
        for record in row.get('points', []):
            identity = (record['path'], record['sha256'])
            if identity in seen:
                continue
            seen.add(identity)
            point, x, native = base.admit_point(world, source, row['metal'], record)
            key = coordinate_hash(x)
            if key not in by_coordinate:
                continue
            candidate = by_coordinate[key]
            require(np.array_equal(np.asarray(candidate['_x']), x), 'Trial/candidate hash collision')
            cell_key = (candidate['candidate_id'], row['metal'])
            if cell_key in cells:
                require(cells[cell_key]['native_energies'] == native, 'Conflicting saved trial components')
            else:
                cells[cell_key] = score_row(helper, candidate['candidate_id'], row['source_id'], row['metal'], record, point, native)
                cells[cell_key].update(reused=True, trial_added_as_candidate=False)
    require(len(candidates) <= 22, 'Four-terminal geometry ceiling exceeded')
    return admitted, by_coordinate


def failed_cache(base, world, slots, by_coordinate, cells):
    blocked = {}
    for row in slots:
        sid, metal = row['source_id'], row['metal']
        source = world['sources'][sid]
        for request in row['requests']:
            require(read(request['record']) == {k: v for k, v in request.items() if k != 'record'},
                    'Saved request receipt differs')
            require(request['source_id'] == sid and request['metal'] == metal
                    and request['recipe'] == base.recipe(world, source, metal), 'Request state/recipe differs')
            x = world['pool'].check_coordinates(world['kinworld'], source, metal,
                                               request['xyz'], request['physical_source'], request['positions'])
            require(request['cache_key'] == base.key_for(world, source, metal, x), 'Request cache identity differs')
            if request['status'] == 'complete':
                require(request['point'] in row['points'], 'Complete request lacks its admitted point record')
                continue
            require(request['status'] == 'unavailable', 'Unknown request outcome')
            key = coordinate_hash(x)
            if key in by_coordinate:
                cid = by_coordinate[key]['candidate_id']
                blocked[(cid, metal)] = {'candidate_id': cid, 'metal': metal, 'request': request['record'],
                                        'reason': 'An earlier exact-state molecular attempt failed; no retry'}
    return {k: value for k, value in blocked.items() if k not in cells}


def summarize(helper, prior, candidates, cells, slots):
    """Keep all old menu values; only evaluate the declared expanded menus."""
    old_menus = prior['comparison']['menus']
    require(len(old_menus) == 15 and {(r['source_scope'], r['menu']) for r in old_menus}
            == {(s, n) for s in (*SOURCES, 'global')
                for n in ('origin', 'prior', 'control', 'candidate', 'union')}, 'Old15 denominator differs')
    expanded = []
    for scope in (*SOURCES, 'global'):
        previous = next(r for r in old_menus if r['source_scope'] == scope and r['menu'] == 'union')
        required = [r for r in slots if scope == 'global' or r['source_id'] == scope]
        ids = list(previous['candidate_ids'])
        for row in required:
            if row['status'] == 'admitted' and row['candidate_id'] not in ids:
                ids.append(row['candidate_id'])
        missing = [{'candidate_id': cid, 'metal': metal} for cid in ids for metal in METALS
                   if cells.get((cid, metal), {}).get('status') != 'complete']
        menu = {'source_scope': scope, 'menu': 'matched_response_union', 'status': 'unavailable',
                'candidate_ids': ids, 'declared_search_slots': len(required),
                'declared_states': 2 * len(ids), 'missing_cells': missing}
        expanded.append(menu)
        if any(r['status'] != 'admitted' for r in required) or missing:
            menu['reason'] = 'A declared terminal search or required paired score is unavailable'
            continue
        selected = {m: min((cells[(cid, m)] for cid in ids),
                           key=lambda row: row['total_energy_kcal_mol']) for m in METALS}
        values = {m: selected[m]['native_energies'] for m in METALS}
        diff = helper.difference(values['Dy'], values['La'])
        old_diff = previous['LBT_Dy_minus_La_components_kcal_mol']
        menu.update(status='available', selected={m: selected[m]['candidate_id'] for m in METALS},
            exact_ties={m: [cid for cid in ids if cells[(cid, m)]['total_energy_kcal_mol']
                            == selected[m]['total_energy_kcal_mol']] for m in METALS},
            selected_native_energies=values, LBT_Dy_minus_La_components_kcal_mol=diff,
            LBT_Dy_minus_La_kcal_mol=math.fsum(diff.values()),
            contrast_change_from_old_union_kcal_mol=math.fsum(diff[k] - old_diff[k] for k in helper.NAMES),
            selected_work_from_old_union_components_kcal_mol={
                m: helper.difference(values[m], previous['selected_native_energies'][m]) for m in METALS})
    local = [row for row in expanded if row['source_scope'] != 'global']
    spread = {'status': 'unavailable'}
    if all(row['status'] == 'available' for row in local):
        signed = math.fsum([local[1]['LBT_Dy_minus_La_kcal_mol'], -local[0]['LBT_Dy_minus_La_kcal_mol']])
        spread = {'status': 'available', '7CCN_minus_7CCO_kcal_mol': signed, 'absolute_spread_kcal_mol': abs(signed)}
    return {'old_comparison_unchanged': prior['comparison'], 'expanded_menus': expanded,
            'expanded_source_spread': spread,
            'contrast_convention': 'Dy minus La; balanced exchange integration with saved Mex12731 is separate',
            'limits': 'Consumed finite source menus, not populations; conditional chemistry, nonstationary endpoints, no affinity claim.'}


def run(cfg, out, result):
    base, world, prior, candidates, cells, slots, parents = inputs(cfg, out)
    result.update(matched_pool_result=cfg['matched_pool_result'], scout_result=cfg['scout_result'],
                  remaining_result=cfg['remaining_result'], old_rows_unchanged=prior['rows'],
                  old_comparison_unchanged=prior['comparison'], old_candidate_count=18, old_cell_count=36,
                  actual_cost={k: 0 for k in base.COST_KEYS})
    admitted, by_coordinate = assemble(base, world, candidates, cells, slots)
    blocked = failed_cache(base, world, slots, by_coordinate, cells)
    missing = [(c['candidate_id'], m) for c in candidates for m in METALS
               if (c['candidate_id'], m) not in cells and (c['candidate_id'], m) not in blocked]
    require(len(missing) <= 4 and all(int(cid.rsplit('_', 1)[1]) >= 18 for cid, _ in missing),
            'Only at most four new terminal opposite-state cells are authorized')
    if missing:
        require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'),
                'New cross-state MACE calls require an allocated GPU')
    public_candidates = [{k: v for k, v in c.items() if k != '_x'} for c in candidates]
    result.update(terminal_slots=admitted, candidates=public_candidates,
                  aliases={alias: c['candidate_id'] for c in candidates for alias in c['aliases']},
                  required_cross_cells=[{'candidate_id': cid, 'metal': m} for cid, m in missing],
                  unavailable_cached_cells=list(blocked.values()), declared_geometries=len(candidates),
                  declared_cells=2 * len(candidates))
    write(out / 'PREFLIGHT.json', {'status': 'passed', 'molecular_calls': 0, 'model_loaded': False,
          'matched_pool_result': cfg['matched_pool_result'], 'declared_search_slots': 4,
          'terminal_slots': admitted, 'old_geometries': 18, 'old_cells': 36,
          'candidates': public_candidates, 'required_cross_cells': result['required_cross_cells'],
          'unavailable_cached_cells': list(blocked.values()), 'cross_call_ceiling': 4})
    result['preflight'] = pin(out / 'PREFLIGHT.json')
    native_model = None
    for cid, metal in missing:
        candidate = next(c for c in candidates if c['candidate_id'] == cid)
        sid = candidate['source_id']
        directory = out / (cid + '_' + metal)
        directory.mkdir()
        evaluator = base.MatchedEvaluator(world, world['sources'][sid], metal, directory, native=native_model)
        row = {'candidate_id': cid, 'source_id': sid, 'metal': metal, 'status': 'unavailable', 'reused': False}
        result['cross_rows'].append(row)
        try:
            require(evaluator.key(candidate['_x']) not in evaluator.cache, 'Unexpected hidden pool cache')
            result['new_cross_cells_attempted'] += 1
            require(result['new_cross_cells_attempted'] <= 4, 'Cross-call ceiling exceeded')
            point, record, hit = evaluator.evaluate(candidate['_x'], compute_force=False)
            require(not hit, 'New cross request unexpectedly returned a cache hit')
            point, x, values = base.admit_point(world, world['sources'][sid], metal, record)
            require(coordinate_hash(x) == candidate['coordinate_sha256'], 'Cross coordinate identity differs')
            raw = read(point['MACE'])
            require(raw['forces_eV_A'] is None and raw['force_convention'] is None,
                    'Pool requested or received an unexpected MACE force')
            row.update(score_row(world['h'], cid, sid, metal, record, point, values))
            cells[(cid, metal)] = row
        except Exception as exc:
            row.update(failure=str(exc), traceback=traceback.format_exc())
        finally:
            native_model = evaluator.native
            row['actual_cost'] = dict(evaluator.cost)
            row['requests'] = evaluator.requests
            for key in base.COST_KEYS:
                result['actual_cost'][key] += evaluator.cost[key]
    require(result['actual_cost']['new_composite_evaluations'] <= 4
            and result['actual_cost']['MACE_calls_attempted'] <= 4
            and result['actual_cost']['ALPB_calls_attempted'] <= 4, 'Molecular call ceiling exceeded')
    result['comparison'] = summarize(world['h'], prior, candidates, cells, admitted)
    result['paired_matrix'] = [{'candidate_id': c['candidate_id'], 'aliases': c['aliases'],
        'cells': {m: cells.get((c['candidate_id'], m), {'status': 'unavailable'}) for m in METALS}} for c in candidates]
    result['status'] = 'complete' if all(r['status'] == 'available' for r in result['comparison']['expanded_menus']) else 'unavailable'
    result['search_costs_unchanged'] = [parent['actual_cost'] for parent in parents]
    result['incremental_call_counts_including_searches'] = {
        k: result['actual_cost'][k] + sum(parent['actual_cost'][k] for parent in parents)
        for k in ('new_composite_evaluations', 'MACE_calls_attempted', 'ALPB_calls_attempted')}
    require(all(v <= 100 for v in result['incremental_call_counts_including_searches'].values()),
            'Declared search-plus-pool ceiling exceeded')
    if native_model is not None:
        result.update(native_model_calls=native_model.calls, native_device=native_model.device)
        require(native_model.calls == result['actual_cost']['MACE_calls_attempted'], 'Warm model call accounting differs')
    require(result['old_rows_unchanged'] == prior['rows']
            and all(cells[(r['candidate_id'], r['metal'])] == r for r in prior['rows']), 'Old36 cells changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit()
            and int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) >= 1, 'Slurm allocation required')
    cfg = json.loads(args.config.read_text())
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick = time.monotonic()
    result = {'protocol_id': PROTOCOL, 'stage': 'pool', 'status': 'running', 'implementation': pin(__file__),
              'config': pin(args.config), 'plan': cfg['plan'], 'slurm_job_id': os.environ['SLURM_JOB_ID'],
              'host': socket.gethostname(), 'declared_search_slots': 4, 'new_cross_cells_attempted': 0,
              'new_MACE_gradients_requested': 0, 'new_GFN2_calls': 0, 'new_MM_calls': 0,
              'new_DFT_calls': 0, 'optimizer_steps': 0, 'cross_rows': []}
    try:
        run(cfg, out, result)
        require(pin(args.config) == result['config'], 'Pool config changed during execution')
    except Exception as exc:
        result.update(status='execution_failure', failure=str(exc), traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic() - tick
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds'] * int(os.environ['SLURM_CPUS_PER_TASK'])
    result['application_allocated_GPU_seconds'] = (result['elapsed_seconds']
        if os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles') else 0.)
    result['peak_host_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['allocation'] = {k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_CPU', 'CUDA_VISIBLE_DEVICES')}
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}), flush=True)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
