#!/usr/bin/env python3
"""Add only admitted range-continuation terminals to the complete old LBT3 pool.

CLI: LBT3_RANGE_POOL.py run --config CONFIG --output FRESH
Config pins: implementation, plan, old_pool_result, range_scout_result,
range_remaining_result. One warm native MACE model, two single-thread GFN2
workers, at most four new opposite-state composite cells. Slurm only.
Old candidates/cells remain immutable; source-local and global reports are
assembled in this same run. No optimization or new gradient requests.
"""
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

POOL_SHA = '505ae7e59ba6abcd2575536f40291bb1c7932a911e21e3aac6268d865885ae54'
RANGE_PROTOCOL = 'lbt3_whole_composite_range04_continuation_v1'
PROTOCOL = 'lbt3_whole_range04_terminal_common_pool_v1'
SOURCES, METALS = ('7CCO', '7CCN'), ('La', 'Dy')


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
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False); handle.write('\n')


def module(name, record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed helper')
    spec = importlib.util.spec_from_file_location(name, record['path'])
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def inputs(config_path):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    cfg = json.loads(config_path.read_text())
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Range pool adapter differs')
    old = read(cfg['old_pool_result']); old_cfg = read(old['config'])
    require(old['status'] == 'complete' and old['stage'] == 'run' and old['slurm_job_id'] == '12541'
            and old['implementation'] == old_cfg['implementation']
            and old['implementation']['sha256'] == POOL_SHA, 'Accepted old pool lineage differs')
    pool = module('lbt3_range_pool_frozen_pool', old['implementation'])
    _, world = pool.inputs(pool.verify(old['config']))
    pool.verify(cfg['plan'])
    candidates, aliases, cells, old_slots, missing = pool.assemble(world)
    require(len(candidates) == 6 and all(s['status'] == 'admitted' for s in old_slots)
            and old['searches'] == old_slots and old['aliases'] == aliases
            and old['candidates'] == [{k: v for k, v in c.items() if k != 'x'} for c in candidates]
            and old['required_cross_cells'] == [{'candidate_id': i, 'metal': m} for i, m in missing],
            'Actual old six-candidate pool replay differs')
    require(len(old['cross_rows']) == len(missing) == 4 and
            {(r['candidate_id'], r['metal']) for r in old['cross_rows']} == set(missing), 'Old cross-state matrix incomplete')
    for row in old['cross_rows']:
        require(row['status'] == 'accepted' and read(row['cell']) == {k: v for k, v in row.items() if k != 'cell'},
                'Old cross-cell receipt differs')
        cells[(row['candidate_id'], row['metal'])] = row
    require(len(cells) == 12 and pool.summarize(candidates, aliases, cells, old_slots) == old['summary'],
            'Old twelve-cell scalar report does not replay exactly')
    # Preserve exact old rows; reuse status is stated at the pool level, not rewritten in old evidence.
    scout, remaining = read(cfg['range_scout_result']), read(cfg['range_remaining_result'])
    require(scout['stage'] == 'scout' and remaining['stage'] == 'remaining'
            and remaining['range_scout_result'] == cfg['range_scout_result'], 'Range stage lineage differs')
    for parent in (scout, remaining):
        parent_cfg = read(parent['config']); pool.verify(parent['implementation'])
        require(parent['status'] in ('complete', 'unavailable') and parent['preflight']['status'] == 'passed'
                and parent['preflight']['declared_states'] == 4 and parent['preflight']['molecular_calls'] == 0
                and parent['preflight']['model_loaded'] is False and read(parent['preflight_result']) == parent['preflight'],
                'Range parent execution/preflight failure cannot admit terminal rows')
        require((parent['status'] == 'complete') == all(r['status'] == 'complete' and r['candidate_admitted'] for r in parent['rows']),
                'Range parent completion status differs from declared search availability')
        require(parent['old_pool_result'] == parent_cfg['old_pool_result'] == cfg['old_pool_result']
                and parent['plan'] == parent_cfg['plan'] == cfg['plan']
                and parent['implementation'] == parent_cfg['implementation']
                and parent['protocol_id'] == RANGE_PROTOCOL
                and parent['profile_result'] == old_cfg['profile_result'] and parent['fd_result'] == old_cfg['fd_result'],
                'Range parent model/source/plan differs')
    require(scout['implementation']['sha256'] == remaining['implementation']['sha256']
            and scout['settings'] == remaining['settings'] == dict(world['parents']['scout']['settings'], all_atom_displacement_A=.40),
            'Only original-source displacement bound may change')
    range_module = module('lbt3_range_pool_guard', scout['implementation'])
    base = module('lbt3_range_pool_original_guard', scout['frozen_scout_implementation'])
    require(scout['frozen_scout_implementation'] == remaining['frozen_scout_implementation'] == world['parents']['scout']['implementation'],
            'Original physical guard lineage differs')
    global np
    import numpy as np
    base.np = range_module.np = np
    require(pin(sys.executable)['sha256'] == world['mcfg']['python']['sha256'], 'Native Python differs')
    slots = [dict(row, parent=cfg[key]) for key, parent in (('range_scout_result', scout), ('range_remaining_result', remaining)) for row in parent['rows']]
    require(len(slots) == 4 and {(s['source_id'], s['metal']) for s in slots} == {(s, m) for s in SOURCES for m in METALS},
            'Four continuation slots required')
    return cfg, pool, world, old, candidates, aliases, cells, old_slots, slots, range_module, base, (scout, remaining)


def add_terminals(pool, world, candidates, aliases, cells, old_slots, slots, range_module, base):
    from scaffold_accommodation import Geometry, SETTINGS as GEOMETRY_SETTINGS
    old_by = {(s['source_id'], s['metal']): s for s in old_slots}
    by_coord = {c['coordinate_sha256']: c for c in candidates}; admitted = []
    active = world['profile']['active_candidate_basis']['mode_indices_selection_order']
    for row in sorted(slots, key=lambda r: (SOURCES.index(r['source_id']), METALS.index(r['metal']))):
        source_id, metal = row['source_id'], row['metal']; source = world['sources'][source_id]
        previous = old_by[(source_id, metal)]; old_search = read(previous['search'])
        result = {'source_id': source_id, 'metal': metal, 'parent': row['parent'], 'status': 'unavailable'}
        if row.get('search') is None:
            result['failure'] = 'No continuation search receipt'; admitted.append(result); continue
        search = read(row['search']); result['search'] = row['search']
        require(search['source_id'] == source_id and search['metal'] == metal, 'Continuation identity differs')
        if row['status'] != 'complete' or not search['candidate_admitted']:
            result['failure'] = 'Continuation unavailable; earlier valid points remain archived, not new candidates'
            admitted.append(result); continue
        require(search['protocol_id'] == RANGE_PROTOCOL and search['status'] == 'finite_proposal_available'
                and search['candidate'] in row['points'] and search['candidate'] == search['trace'][-1]['point']
                and search['continuation_start'] == search['trace'][0]['point'] == previous['point']
                and search['previous_search'] == previous['search'] and search['source_rebased'] is False
                and search['coordinate_scaling_T'] == old_search['coordinate_scaling_T']
                and search['authoritative_origin'] == old_search['authoritative_origin']
                and search['active_indices'] == active and search['active_names'] == base.NAMES
                and search['initial_history_length'] == search['initial_new_evaluations'] == 0
                and len(search['trace']) == search['accepted_steps'] + 1
                and search['origin_guard']['pass'] and search['continuation_start_guard']['pass']
                and search['final_geometry']['pass'] and all(t['geometry']['pass'] for t in search['trace']),
                'Original-reference continuation/admission binding differs')
        point = read(search['candidate']); x, components = pool.point_components(world, source, metal, point)
        q = np.asarray(search['full_q']); replay, _, _, _ = source['kin'].evaluate(q)
        require(q.shape == (64,) and np.isfinite(q).all()
                and np.count_nonzero(q[[k for k in range(64) if k not in active]]) == 0
                and np.allclose(replay, x, rtol=0, atol=1e-12), 'New terminal kinematic scope/replay differs')
        guard = range_module.range_guard(base, source, active, Geometry, GEOMETRY_SETTINGS, metal)
        check = guard.check(q, x, x)
        require(check['pass'] and check['source_baseline'] == source['prepared']['physical_atoms']
                and check['source_displacement_bound_A'] == .40, 'New terminal fails original-source range guards')
        key = hashlib.sha256(np.where(x == 0., 0., x).astype('<f8').tobytes()).hexdigest()
        if key not in by_coord:
            candidate = {'candidate_id': 'geometry_%02d' % len(candidates), 'coordinate_sha256': key,
                'source_id': source_id, 'physical_source': point['physical_source'], 'aliases': [], 'x': x}
            candidates.append(candidate); by_coord[key] = candidate
        candidate = by_coord[key]; require(np.array_equal(candidate['x'], x), 'Exact coordinate hash collision')
        alias = source_id + ':range_terminal_' + metal
        candidate['aliases'].append(alias); aliases[alias] = candidate['candidate_id']
        cell_key = (candidate['candidate_id'], metal)
        if cell_key in cells:
            require(cells[cell_key]['components'] == components, 'Conflicting exact-state reuse; old cell cannot be replaced')
        else:
            cells[cell_key] = {'candidate_id': candidate['candidate_id'], 'metal': metal, 'status': 'accepted',
                'components': components, 'reused': True, 'evidence': [{'terminal_point': search['candidate'], 'search': row['search'], 'parent': row['parent']}]}
        work = pool.difference(components, world['origin_components'][(source_id, metal)])
        increment = pool.difference(components, read(previous['point'])['components'])
        require(work == search['work_from_origin'] and increment == search['work_from_continuation_start']
                and increment['composite_kcal_mol'] <= 1e-8, 'Continuation work/reference differs')
        result.update(status='admitted', candidate_id=candidate['candidate_id'], point=search['candidate'],
            previous_point=previous['point'], physical_check=check, own_work_from_original=work,
            own_work_from_continuation_start=increment, termination=search['termination'], accepted_steps=search['accepted_steps'],
            final_force_diagnostic=search['final_force_diagnostic'], boundary_flags=search['preserved_boundary_flags'])
        admitted.append(result)
    missing = [(c['candidate_id'], m) for c in candidates for m in METALS if (c['candidate_id'], m) not in cells]
    require(len(candidates) <= 10 and len(missing) <= 4, 'Expanded candidate/cross-cell ceiling exceeded')
    return admitted, missing


def summarize(pool, old, candidates, aliases, cells, slots):
    result = {'prior_summary': old['summary'], 'source_local_menus': {}, 'source_spread': {},
        'global_menu_caveat': 'The shared global menu removes source choice by construction; only source-local contrasts diagnose remaining source sensitivity.'}
    complete = all(s['status'] == 'admitted' for s in slots)
    for source in SOURCES:
        prior = old['summary']['source_local_pools'][source]['plus_both_terminals']
        names = [source + ':origin', source + ':terminal_La', source + ':terminal_Dy',
                 source + ':range_terminal_La', source + ':range_terminal_Dy']
        expanded = pool.pool_summary([aliases[n] for n in names], cells) if all(n in aliases for n in names) else {
            'status': 'unavailable', 'reason': 'A declared continuation terminal is unavailable'}
        local = {'prior': prior, 'expanded': expanded}
        if expanded['status'] == 'available':
            local['contrast_change_kcal_mol'] = expanded['independently_selected_Dy_minus_La']['composite_kcal_mol'] - prior['independently_selected_Dy_minus_La']['composite_kcal_mol']
            local['own_metal_work'] = {m: pool.difference(cells[(expanded['selected'][m], m)]['components'], cells[(prior['selected'][m], m)]['components']) for m in METALS}
        result['source_local_menus'][source] = local
    for menu in ('prior', 'expanded'):
        pair = [result['source_local_menus'][s][menu] for s in SOURCES]
        if all(p['status'] == 'available' for p in pair):
            a, b = [p['independently_selected_Dy_minus_La']['composite_kcal_mol'] for p in pair]
            result['source_spread'][menu] = {'7CCN_minus_7CCO_kcal_mol': b-a, 'absolute_spread_kcal_mol': abs(b-a)}
        else:
            result['source_spread'][menu] = {'status': 'unavailable'}
    if 'absolute_spread_kcal_mol' in result['source_spread']['expanded']:
        result['source_spread_change_kcal_mol'] = result['source_spread']['expanded']['absolute_spread_kcal_mol'] - result['source_spread']['prior']['absolute_spread_kcal_mol']
    prior = old['summary']['expanded_global_common_pool']
    expanded = pool.pool_summary([c['candidate_id'] for c in candidates], cells) if complete else {
        'status': 'unavailable', 'reason': 'Four admitted continuation searches are required; prior pool remains available separately'}
    result.update(prior_global_pool=prior, expanded_global_pool=expanded)
    if expanded['status'] == 'available':
        result['global_change'] = {'contrast_change_kcal_mol': expanded['independently_selected_Dy_minus_La']['composite_kcal_mol'] - prior['independently_selected_Dy_minus_La']['composite_kcal_mol'],
            'selection_changed': {m: expanded['selected'][m] != prior['selected'][m] for m in METALS},
            'own_metal_work': {m: pool.difference(cells[(expanded['selected'][m], m)]['components'], cells[(prior['selected'][m], m)]['components']) for m in METALS}}
    result['paired_matrix'] = [{'candidate_id': c['candidate_id'], 'aliases': c['aliases'],
        'cells': {m: cells.get((c['candidate_id'], m), {'status': 'unavailable'}) for m in METALS}} for c in candidates]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('stage', choices=['run'])
    for key in ('config', 'output'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    cfg, pool, world, old, candidates, aliases, cells, old_slots, slots, range_module, base, parents = inputs(args.config.resolve())
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    result = {'protocol_id': PROTOCOL, 'stage': 'run', 'status': 'unavailable', 'implementation': pin(__file__), 'config': pin(args.config),
        'plan': cfg['plan'], 'old_pool_result': cfg['old_pool_result'], 'range_scout_result': cfg['range_scout_result'],
        'range_remaining_result': cfg['range_remaining_result'], 'pool_helper': old['implementation'],
        'new_gradients_requested': 0, 'optimizer_steps': 0, 'prior_candidate_count': 6, 'prior_cell_count': 12,
        'declared_continuations': 4, 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
        'actual_cost': {'new_composite_evaluations': 0, 'MACE_calls_attempted': 0, 'GFN2_calls_attempted': 0,
                        'model_load_seconds': 0., 'MACE_call_seconds': 0., 'GFN2_call_seconds': 0.}, 'cross_rows': [],
        'prior_summary': old['summary'], 'interpretation_limits': 'Conditional intact-peptide energy response; no affinity, equilibrium, classifier or Mex composite inference.'}
    write(out / 'START.json', result)
    try:
        admitted, missing = add_terminals(pool, world, candidates, aliases, cells, old_slots, slots, range_module, base)
        result.update(searches=admitted, aliases=aliases, candidates=[{k: v for k, v in c.items() if k != 'x'} for c in candidates],
            required_cross_cells=[{'candidate_id': i, 'metal': m} for i, m in missing])
        write(out / 'POOL.json', {k: result[k] for k in ('searches', 'aliases', 'candidates', 'required_cross_cells')})
        result['pool'] = pin(out / 'POOL.json')
        rows = pool.execute(world, candidates, missing, out, result) if all(s['status'] == 'admitted' for s in admitted) else []
        for row in rows:
            require(read(row['cell']) == {k: v for k, v in row.items() if k != 'cell'}, 'New cross-cell receipt differs')
            if row['status'] == 'accepted':
                cells[(row['candidate_id'], row['metal'])] = row
        result['cross_rows'] = rows; result['summary'] = summarize(pool, old, candidates, aliases, cells, admitted)
        result['status'] = 'complete' if result['summary']['expanded_global_pool']['status'] == 'available' else 'unavailable'
        result['continuation_costs'] = [p['actual_cost'] for p in parents]
        result['incremental_call_counts_including_searches'] = {key: result['actual_cost'][key] + sum(p['actual_cost'][key] for p in parents)
            for key in ('new_composite_evaluations', 'MACE_calls_attempted', 'GFN2_calls_attempted')}
        require(result['incremental_call_counts_including_searches']['MACE_calls_attempted'] <= 100 and
                result['incremental_call_counts_including_searches']['GFN2_calls_attempted'] <= 200, 'Declared whole-extension call ceiling exceeded')
        pool.verify(result['config'])
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic() - started
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds'] * int(os.environ['SLURM_CPUS_PER_TASK'])
    result['application_allocated_GPU_seconds'] = result['elapsed_seconds']
    result['peak_host_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['allocation'] = {k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_CPU', 'CUDA_VISIBLE_DEVICES')}
    write(out / 'RESULT.json', result); print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}), flush=True)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
