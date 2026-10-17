#!/usr/bin/env python3
"""Mex common5 to common7, same Cartesian composite, at most two cross cells."""
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

PROTOCOL = 'mex_full_cartesian_common7_v1'
METALS = ('La', 'Dy')
OLD_ALIASES = ['origin', 'La_pilot', 'Dy_pilot', 'La_terminal', 'Dy_terminal']
ALIASES = OLD_ALIASES+['La_cartesian_terminal', 'Dy_cartesian_terminal']


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


def coordinate_hash(ev, x):
    return hashlib.sha256(ev.np.asarray(x, dtype='<f8').tobytes(order='C')).hexdigest()


def update_work(ev, row):
    if any(row['cells'].get(m, {}).get('status') != 'complete' for m in METALS):
        row['status'] = 'unavailable'; return
    components = {}; totals = {}
    for metal in METALS:
        parts = ev.base.relative(ev.old.response, ev.old.old, row['cells'][metal]['native_energies'],
                                  ev.old.response.origin_native(metal, ev.old.world))
        components[metal] = {k: v for k, v in parts.items() if k != 'composite_kcal_mol'}
        totals[metal] = parts['composite_kcal_mol']
    row.update(status='complete', component_work_from_original_source_kcal_mol=components,
        total_work_from_original_source_kcal_mol=totals, common_work_kcal_mol=(totals['La']+totals['Dy'])/2,
        Dy_minus_La_work_kcal_mol=totals['Dy']-totals['La'])


def summarize(oldpool, ev, matrix, aliases, slots):
    missing_aliases = [a for a in aliases if not any(a in row['aliases'] for row in matrix)]
    needed = [row for row in matrix if any(a in aliases for a in row['aliases'])]
    missing = [{'coordinate_sha256': row['coordinate_sha256'], 'metal': m} for row in needed for m in METALS
               if row['cells'].get(m, {}).get('status') != 'complete']
    result = {'status': 'unavailable', 'declared_menu': aliases, 'missing_aliases': missing_aliases,
              'missing_cells': missing, 'declared_search_slots': len(slots)}
    if missing_aliases or missing or any(s['status'] != 'admitted' for s in slots):
        result['reason'] = 'A declared terminal or required paired state is unavailable'; return result
    result.update(oldpool.summary(ev.base, ev.old, matrix, aliases))
    result['exact_ties'] = {m: [a for a in aliases if next(r for r in matrix if a in r['aliases'])[
        'total_work_from_original_source_kcal_mol'][m] == result['selections'][m]['work_from_original_source_kcal_mol']]
        for m in METALS}
    return result


def cached_cell(ev, metal, x):
    record = ev.old.cache[ev.key(metal, x)]; point = read(record)
    require(point['metal'] == metal and ev.np.array_equal(ev.old.coordinates(point), x), 'Cached physical state differs')
    ev.package(point, record, x, True)
    return {'status': 'complete', 'metal': metal, 'point': record, 'native_energies': point['native_energies'],
            'assembled_gradient': point['assembled_gradient'], 'physical': point['physical'], 'reused': True}


def baseline(ev, oldpool, guard, out):
    require(ev.pool['slurm_job_id'] == '12731' and ev.pool['status'] == 'complete'
            and ev.pool['declared_menu'] == OLD_ALIASES, 'Exact completed common5 required')
    matrix = []; by_coordinate = {}
    for saved in ev.pool['matrix']:
        x = ev.old.coordinates(read(saved['cells']['La']['point'])); key = coordinate_hash(ev, x)
        require(key == saved['coordinate_sha256'] and key not in by_coordinate, 'Old exact coordinate identity differs')
        row = copy.deepcopy(saved); row['_x'] = x; matrix.append(row); by_coordinate[key] = row
        for metal in METALS:
            cell = cached_cell(ev, metal, x)
            require(cell['native_energies'] == saved['cells'][metal]['native_energies'], 'Old paired energy differs')
            row['cells'][metal]['status'] = 'complete'
        update_work(ev, row)
        require(all(row[k] == saved[k] for k in ('component_work_from_original_source_kcal_mol',
                'total_work_from_original_source_kcal_mol', 'common_work_kcal_mol', 'Dy_minus_La_work_kcal_mol')),
                'Old work replay differs')
        check = guard.check(x); require(check['pass'] and check['state_compatible'], 'Old pool physical/state admission differs')
    require(len(matrix) == 5 and sum(len(r['cells']) for r in matrix) == 10, 'Old5/10 denominator differs')
    summary = copy.deepcopy(summarize(oldpool, ev, matrix, OLD_ALIASES, []))
    require(summary['status'] == 'available' and all(summary[k] == v for k, v in ev.pool['common5'].items()),
            'Original common5 selection/work replay differs')
    partial = copy.deepcopy(matrix); partial[0]['cells'].pop('La')
    test = summarize(oldpool, ev, partial, OLD_ALIASES, [])
    require(test['status'] == 'unavailable' and len(test['missing_cells']) == 1, 'Missing actual-cell test failed')
    record = {'status': 'passed', 'new_molecular_calls': 0, 'old_common5': summary, 'old_geometries': 5, 'old_cells': 10,
              'missing_cell_test': {'test_only': True, 'removed_actual_cell': [matrix[0]['coordinate_sha256'], 'La'], 'result': test}}
    write(out/'BASELINE_REPLAY.json', record)
    return matrix, by_coordinate, summary, pin(out/'BASELINE_REPLAY.json')


def import_request(ev, metal, record):
    request = read(record)
    require(request['identity_sha256'] == ev.identity_hash and request['state'] == ev.state_metadata[metal],
            'Saved Cartesian request model/state differs')
    previous = ev.old.cache.get(request['cache_key'])
    if previous is not None:
        require(read(previous)['native_energies'] == request['native_energies'],
                'Conflicting exact cached native components')
    return ev.restore_cartesian(record)


def restore(search_module, ev, metal, provenance):
    request = read(provenance['cartesian_request'])
    require(provenance['point'] == request['point'] and provenance['cache_key'] == request['cache_key']
            and provenance['identity_sha256'] == request['identity_sha256'] == ev.identity_hash
            and provenance['state'] == request['state'] == ev.state_metadata[metal], 'Cartesian request/provenance differs')
    import_request(ev, metal, provenance['cartesian_request'])
    return search_module.SearchAdapter(ev, metal).read_point(provenance)


def admit_searches(cfg, ev, search_module, optimizer, guard, matrix, by_coordinate):
    require(len(cfg['search_results']) == 2, 'Two ordered La/Dy search result pins required')
    slots = []; failures = []
    for metal, parent_pin in zip(METALS, cfg['search_results']):
        parent = read(parent_pin); scfg = read(parent['config'])
        verify(parent['implementation']); verify(scfg['implementation'])
        require(parent['implementation']['sha256'] == scfg['implementation']['sha256'] == cfg['search_implementation']['sha256']
                and parent['stage'] == 'run' and parent['metal'] == metal, 'Search source/state differs')
        require(read(scfg['evaluator_config']) == ev.cfg and scfg['optimizer']['sha256'] == cfg['optimizer']['sha256']
                and scfg['protocol'] == search_module.PROTOCOL and scfg['scope'] == search_module.SCOPE,
                'Search model/optimizer/scope differs')
        for name in ('optimizer', 'plan', 'specification', 'evaluator'):
            verify(scfg[name])
        search_module.gate(ev, scfg)
        slot = {'metal': metal, 'status': 'unavailable', 'search_result': parent_pin}; slots.append(slot)
        # Import complete rejected and accepted trials with explicit Cartesian
        # admission. No new point is passed through the old q7 map.
        for request_pin in parent.get('requests', []):
            request = read(request_pin)
            require(request['state'] == ev.state_metadata[metal], 'Saved request target state differs')
            import_request(ev, metal, request_pin)
        for request in parent.get('native_requests', []):
            point = read(request['point'])
            require(point['metal'] == request['metal'] == metal, 'Native request state differs')
            if point['status'] == 'complete':
                continue
            require(point['status'] == 'unavailable', 'Unknown native failure state')
            item = {'metal': metal, 'point': request['point'], 'reason': 'Earlier molecular attempt failed; no retry'}
            if 'physical' in point:
                physical = read(point['physical'])
                require([a['id'] for a in physical] == ev.ids, 'Failed request physical map differs')
                item['coordinate_sha256'] = coordinate_hash(ev, [a['xyz_A'] for a in physical])
            failures.append(item)
        if 'search' not in parent:
            slot['reason'] = 'Search result lacks a terminal receipt'; continue
        search = read(parent['search']); slot.update(search=parent['search'], termination=search['termination'],
            full_force=search['final_full_force'], accepted_steps=search['accepted_steps'],
            distinct_evaluations=search['distinct_evaluations'], converged=search['converged'],
            work_from_start=search['work_from_start'], final_geometry=search['final_geometry'])
        checkpoint = read(search['checkpoint'])
        require(checkpoint['point'] == search['terminal_point'] and checkpoint['initial_point'] == search['start_point']
                and checkpoint['accepted_steps'] == search['accepted_steps']
                and checkpoint['distinct_evaluations'] == search['distinct_evaluations']
                and checkpoint['settings'] == search['settings'] and checkpoint['limits'] == search['limits'],
                'Final cumulative checkpoint differs')
        for provenance in checkpoint['evaluated_points']:
            restore(search_module, ev, metal, provenance)
        for request_pin in checkpoint['evaluation_requests']:
            import_request(ev, metal, request_pin)
        require(search['settings'] == optimizer.SETTINGS and search['limits'] == search_module.SCOPE['search_limits']
                and 0 <= search['accepted_steps'] <= 80 and 0 <= search['distinct_evaluations'] <= 200,
                'Search rule/budget differs')
        initial = restore(search_module, ev, metal, search['start_point'])
        require(ev.np.array_equal(initial['positions_A'], ev.starts[metal])
                and initial['native_energies'] == read(ev.cfg['starts'][metal])['native_energies'],
                'Search did not start at exact original restricted terminal')
        if parent['status'] != 'complete' or search['status'] != 'complete' or not search['physical_state_admitted']:
            slot['reason'] = 'Execution or bound-state admission failed'; continue
        terminal = restore(search_module, ev, metal, search['terminal_point'])
        x = ev.np.asarray(terminal['positions_A']); check = guard.check(x)
        require(check == search['final_geometry'] and check['pass'] and check['state_compatible']
                and optimizer.force_metrics(terminal['gradient_kcal_mol_A']) == search['final_full_force'],
                'Terminal physical/full-force replay differs')
        key = coordinate_hash(ev, x); row = by_coordinate.get(key)
        if row is None:
            row = {'aliases': [], 'coordinate_sha256': key, 'cells': {}, 'status': 'unavailable', '_x': x}
            matrix.append(row); by_coordinate[key] = row
        require(ev.np.array_equal(row['_x'], x), 'Physical coordinate hash collision')
        alias = metal+'_cartesian_terminal'; require(alias not in row['aliases'], 'Duplicate new alias')
        row['aliases'].append(alias); cell = cached_cell(ev, metal, x)
        if metal in row['cells']:
            require(row['cells'][metal]['native_energies'] == cell['native_energies'], 'Exact own-state energy conflict')
        else:
            row['cells'][metal] = cell
        slot.update(status='admitted', alias=alias, coordinate_sha256=key, terminal_point=search['terminal_point'])
    require(len(matrix) <= 7, 'Two-new-geometry ceiling exceeded')
    for row in matrix:
        for metal in METALS:
            if metal not in row['cells'] and ev.key(metal, row['_x']) in ev.old.cache:
                row['cells'][metal] = cached_cell(ev, metal, row['_x'])
        update_work(ev, row)
    return slots, failures


def run(cfg, out, result, stage):
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Pool source differs')
    verify(cfg['plan']); verify(cfg['specification'])
    ecfg = read(cfg['evaluator_config']); verify(cfg['evaluator_implementation']); verify(ecfg['implementation'])
    require(ecfg['implementation']['sha256'] == cfg['evaluator_implementation']['sha256'], 'Evaluator source differs')
    em = module('mex_pool_evaluator', ecfg['implementation'])
    search_module = module('mex_pool_search', cfg['search_implementation'])
    optimizer = module('mex_pool_optimizer', cfg['optimizer'])
    ev = em.MexEvaluator(ecfg, out/'evaluations'); result['_ev'] = ev
    guard, _ = search_module.build_guard(ev, optimizer, out)
    oldpool = module('mex_pool_old_common5', ev.cfg['pool_implementation'])
    matrix, by_coordinate, old, replay = baseline(ev, oldpool, guard, out)
    result.update(old_common5=old, old_pool_result=ev.cfg['pool_result'], baseline_replay=replay,
                  identity_sha256=ev.identity_hash, physical_source_ids=ev.ids)
    if stage == 'baseline':
        result.update(status='complete', model_initialized=False, declared_geometries=5, declared_cells=10); return
    slots, failures = admit_searches(cfg, ev, search_module, optimizer, guard, matrix, by_coordinate)
    result.update(terminal_slots=slots, failed_exact_requests=failures, declared_geometries=len(matrix), declared_cells=2*len(matrix))
    all_admitted = all(s['status'] == 'admitted' for s in slots)
    for row in matrix:
        for metal in METALS:
            if metal in row['cells']:
                continue
            cell = {'status': 'unavailable', 'metal': metal, 'coordinate_sha256': row['coordinate_sha256']}
            result['cross_rows'].append(cell)
            if not all_admitted:
                cell['reason'] = 'Both declared search slots must be admitted'; continue
            prior = [f for f in failures if f['metal'] == metal and
                     (f.get('coordinate_sha256') is None or f['coordinate_sha256'] == row['coordinate_sha256'])]
            if prior:
                cell.update(reason='Earlier exact or unlocated failed state; no retry', failed_requests=prior); continue
            require(any(a.endswith('_cartesian_terminal') for a in row['aliases']), 'An old paired score is missing')
            require(result['new_cross_cells_attempted'] < 2, 'Two cross-cell ceiling exceeded')
            require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'), 'Allocated GPU required')
            check = guard.check(row['_x']); require(check['pass'] and check['state_compatible'], 'Cross physical/state gate failed')
            result['new_cross_cells_attempted'] += 1
            try:
                value = ev.evaluate(row['_x'], metal, guard=check)
                saved = read(value['provenance']['point']); admitted = cached_cell(ev, metal, row['_x'])
                require(saved['native_energies'] == admitted['native_energies'], 'Cross native replay differs')
                cell.update(admitted, cartesian_request=ev.requests[-1], reused=False); row['cells'][metal] = cell
            except Exception as exc:
                cell.update(failure=repr(exc), traceback=traceback.format_exc())
        update_work(ev, row)
    # Original aliases select exactly the same cells, despite added alias names.
    original = summarize(oldpool, ev, matrix, OLD_ALIASES, [])
    for metal in METALS:
        require(original['selections'][metal]['point'] == old['selections'][metal]['point'] and
                original['selections'][metal]['work_from_original_source_kcal_mol'] == old['selections'][metal]['work_from_original_source_kcal_mol'],
                'Original selected physical state changed')
    new = summarize(oldpool, ev, matrix, ALIASES, slots)
    result['common7'] = new
    if new['status'] == 'available':
        values = {}
        for metal in METALS:
            a = read(new['selections'][metal]['point'])['native_energies']; b = read(old['selections'][metal]['point'])['native_energies']
            values[metal] = ev.base.relative(ev.old.response, ev.old.old, a, b)
        result['selected_work_from_common5'] = values
        result['common5_to_common7_Dy_minus_La_change_kcal_mol'] = values['Dy']['composite_kcal_mol']-values['La']['composite_kcal_mol']
    result['matrix'] = [{k: v for k, v in row.items() if k != '_x'} for row in matrix]
    result['status'] = 'complete' if all_admitted and new['status'] == 'available' else 'unavailable'
    require(ev.new_evaluations <= 2 and ev.counters['new_MACE_calls_attempted'] <= 2
            and ev.counters['new_GFN2_calls_attempted'] <= 4 and ev.counters['new_MM_configurations_attempted'] <= 2
            and ev.counters['new_MM_group_queries_attempted'] <= 8 and ev.counters['new_ALPB_calls_attempted'] <= 2,
            'Declared composite/component ceiling exceeded')
    lines = ['# Mex full-coordinate common7', '', 'Status: '+result['status']+'.', '',
             'Original common5 remains explicit; both terminal slots and full residual forces are retained.',
             'Conditional EF3 exchange with three Nd spectators; no stationary minimum or affinity claim.', '']
    if new['status'] == 'available':
        lines += ['Additional work from old common5 (kcal/mol):', '',
                  '- La: '+str(result['selected_work_from_common5']['La']['composite_kcal_mol']),
                  '- Dy: '+str(result['selected_work_from_common5']['Dy']['composite_kcal_mol']),
                  '- Dy-minus-La change: '+str(result['common5_to_common7_Dy_minus_La_change_kcal_mol']), '']
    (out/'REPORT.md').write_text('\n'.join(lines)); result['report'] = pin(out/'REPORT.md')


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--stage', choices=('baseline', 'pool'), required=True)
    p.add_argument('--config', type=Path, required=True); p.add_argument('--output', type=Path, required=True); args = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    cfg = json.loads(args.config.read_text()); out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic()
    result = {'protocol_id': PROTOCOL, 'status': 'execution_failure', 'stage': args.stage,
              'implementation': pin(__file__), 'config': pin(args.config), 'plan': cfg['plan'],
              'specification': cfg['specification'], 'search_results': cfg.get('search_results'),
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'new_cross_cells_attempted': 0, 'cross_rows': [],
              'new_DFT_calls': 0, 'new_global_GFN_calls': 0, 'optimizer_steps': 0,
              'declared_terminal_slots': 0 if args.stage == 'baseline' else 2}
    try:
        run(cfg, out, result, args.stage); verify(result['config']); verify(cfg['specification'])
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    ev = result.pop('_ev', None)
    if ev is not None:
        result.update(counters=ev.counters, requests=ev.requests, native_requests=ev.old.requests, native_points=ev.old.points)
    result['elapsed_seconds'] = time.monotonic()-tick; write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
