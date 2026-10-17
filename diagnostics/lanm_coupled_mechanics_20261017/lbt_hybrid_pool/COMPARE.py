#!/usr/bin/env python3
"""Fixed complete Module2 LBT pool under two qualified hybrid partitions."""
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

PROTOCOL = 'lbt_complete_module2_pool_two_hybrid_regions_v1'
QUALIFIER_SHA = 'ca137dd0fa0647268e5be13d21d20b414173767d61f0f519b9a1b693d32cf9df'
POOL_SHA = 'f4663b24190efcbe1707d604dc30a9140b28dfecb048335ec723e961e4401872'
SOURCES = ('7CCO', '7CCN')
REGIONS = ('baseline', 'expanded')
METALS = ('La', 'Dy')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        h = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': h}


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


def same_bytes(a, b):
    verify(a)
    verify(b)
    return a['sha256'] == b['sha256']


def module(name, p):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def qualified_hybrid(cfg, out, result):
    q = read(cfg['qualification_result'])
    qcfg = read(q['config'])
    require(q['status'] == 'complete' and q['stage'] == 'run' and q['qualification_passed']
            and q['implementation']['sha256'] == cfg['qualification_implementation']['sha256'] == QUALIFIER_SHA
            and same_bytes(qcfg['implementation'], cfg['qualification_implementation'])
            and same_bytes(qcfg['plan'], cfg['plan']), 'Complete frozen hybrid qualification required')
    verify(q['implementation'])
    qualifier = module('complete_hybrid_pool_qualification', cfg['qualification_implementation'])
    require(q['protocol'] == qualifier.PROTOCOL and q['declared_new_composites'] == 24
            and q['declared_reused_centers'] == 3 and len(q['checks']) == 3, 'Qualification denominator differs')
    d = out / 'hybrid_parent_replay'
    d.mkdir()
    scratch = {}
    base, ev, centers, frozen = qualifier.load_inputs(qcfg, d, scratch)
    result['_evaluator'] = ev
    require(q['identity_sha256'] == ev.identity_hash and read(q['directions']) == frozen,
            'Qualifier model/directions differ')
    slots = q['declared_displaced_cells']
    keys = [(r['source_id'], r['region'], r['metal'], r['direction'], r['h_A'], r['sign']) for r in slots]
    expected = [(s, r, m, n, h, sign) for s, r, m in qualifier.CASES
                for n in qualifier.FD['directions'] for h in qualifier.FD['steps_A'] for sign in ('plus', 'minus')]
    require(keys == expected and all(r['status'] == 'complete' for r in slots), 'All 24 qualified cells required')
    points = {}
    for key, slot in zip(keys, slots):
        s, r, m, n, h, sign = key
        p = qualifier.restore(ev, slot['point'])
        v = ev.np.asarray(frozen['regions'][r]['directions'][n]['vector'])
        x = ev.origins[(s, r)] + (1 if sign == 'plus' else -1) * h * v
        require((p['source_id'], p['region'], p['metal']) == (s, r, m)
                and ev.np.array_equal(p['positions_A'], x) and p['provenance']['cache_key'] == slot['cache_key'],
                'Qualified cell state/geometry differs')
        points[key] = p
    replay = []
    for s, r, m in qualifier.CASES:
        case = {'source_id': s, 'region': r, 'metal': m, 'status': 'complete',
                'center': centers[(s, r, m)]['provenance']['point'], 'directions': [], 'pass': True}
        for n in qualifier.FD['directions']:
            data = frozen['regions'][r]['directions'][n]
            pairs = [(h, points[(s, r, m, n, h, 'plus')], points[(s, r, m, n, h, 'minus')])
                     for h in qualifier.FD['steps_A']]
            check = qualifier.checks(base, ev, centers[(s, r, m)], ev.np.asarray(data['vector']), pairs)
            require(check['pass'], 'Qualified derivative no longer passes')
            case['directions'].append(dict(check, direction=n, physical_id=data['physical_id']))
        replay.append(case)
    require(replay == q['checks'], 'Actual remaining qualification checks do not replay')
    cost = q['actual_cost']
    require(cost['composites_attempted'] == cost['MACE_calls_attempted'] == cost['MM_configurations_attempted'] == 24
            and cost['GFN_calls_attempted'] == 48 and cost['MM_group_queries_attempted'] == 96
            and cost['ALPB_calls_attempted'] <= 24, 'Qualifier call accounting differs')
    result['hybrid_parent_replay'] = write(out / 'HYBRID_PARENT_REPLAY.json', {
        'status': 'passed', 'scout_and_origins_replay': scratch['preflight'],
        'remaining_qualification': cfg['qualification_result'], 'checks': replay, 'new_molecular_calls': 0})
    result['parent_costs'] = dict(scratch['parent_costs'], remaining_qualification=cost)
    return base, ev


def module2_pool(cfg, out, ev, result):
    parent = read(cfg['module2_pool_result'])
    pcfg = read(parent['config'])
    require(parent['status'] == 'complete' and parent['stage'] == 'run'
            and parent['implementation']['sha256'] == pcfg['implementation']['sha256'] == POOL_SHA
            and same_bytes(pcfg['plan'], cfg['plan']), 'Complete Module2 pool required')
    verify(parent['implementation'])
    pool = module('hybrid_module2_pool_replay', pcfg['implementation'])
    require(parent['protocol_id'] == pool.PROTOCOL and parent['declared_terminal_slots'] == 4, 'Module2 protocol differs')
    cart = pool.module('hybrid_module2_cartesian', pcfg['cartesian_implementation'])
    ccfg = read(pcfg['cartesian_config'])
    require(same_bytes(ccfg['implementation'], pcfg['cartesian_implementation']), 'Module2 Cartesian source differs')
    d = out / 'module2_replay'
    d.mkdir()
    ld = d / 'source_admission'
    ld.mkdir()
    loaded = cart.load_inputs(ccfg, ld)
    pool.np = ev.np
    bd = d / 'old_baseline'
    bd.mkdir()
    candidates, cells, by_coordinate, baselines, baseline_pin = pool.replay_baseline(loaded, bd)
    old_ids = [c['candidate_id'] for c in candidates]
    slots, blocked = pool.admit_runs(pcfg, cart, ccfg, loaded, d, candidates, cells, by_coordinate)
    require(len(slots) == 4 and all(s['status'] == 'admitted' for s in slots) and slots == parent['terminal_slots']
            and not blocked and len(candidates) == parent['declared_geometries'] <= 26,
            'Module2 complete physical terminal admission differs')
    metadata = [{k: v for k, v in c.items() if k != '_x'} for c in candidates]
    require(metadata == parent['candidates'] and len(parent['paired_matrix']) == len(candidates)
            and parent['declared_cells'] == 2 * len(candidates), 'Module2 candidate denominator/order differs')
    base, world, helper = loaded['base'], loaded['world'], loaded['world']['h']
    for candidate, matrix_row in zip(candidates, parent['paired_matrix']):
        cid, sid = candidate['candidate_id'], candidate['source_id']
        require(matrix_row['candidate_id'] == cid and matrix_row['source_id'] == sid
                and matrix_row['aliases'] == candidate['aliases'] and set(matrix_row['cells']) == set(METALS),
                'Module2 matrix identity differs')
        for metal in METALS:
            row = matrix_row['cells'][metal]
            require(row['status'] == 'complete' and (row['candidate_id'], row['source_id'], row['metal']) == (cid, sid, metal),
                    'Required Module2 paired state missing')
            if 'point' in row:
                point, x, native = base.admit_point(world, world['sources'][sid], metal, row['point'], False)
                require(ev.np.array_equal(x, candidate['_x']), 'Module2 POINT coordinates differ')
            else:
                require((cid, metal) in cells, 'Unreplayed old raw MACE cell')
                native = cells[(cid, metal)]['native_energies']
            parts = helper.native_components(native)
            require(native == row['native_energies'] and parts == row['component_energies_kcal_mol']
                    and math.fsum(parts.values()) == row['total_energy_kcal_mol'], 'Module2 component replay differs')
            if (cid, metal) in cells:
                require(cells[(cid, metal)]['native_energies'] == native, 'Conflicting exact Module2 state')
            cells[(cid, metal)] = copy.deepcopy(row)
    require(pool.comparison(helper, candidates, cells, slots, baselines) == parent['comparison'],
            'Module2 selections/components/source spread do not replay')
    require(len(old_ids) == 22 and 22 <= len(candidates) <= 26
            and all(11 <= sum(c['source_id'] == s for c in candidates) <= 13 for s in SOURCES),
            'Complete inherited/candidate source denominator differs')
    result.update(module2_parent=cfg['module2_pool_result'], module2_comparison=parent['comparison'],
                  module2_terminal_slots=slots, parent_costs=dict(result['parent_costs'], module2_pool=parent['actual_cost']))
    result['module2_replay'] = write(out / 'MODULE2_REPLAY.json', {
        'status': 'passed', 'parent': cfg['module2_pool_result'], 'baseline_replay': baseline_pin,
        'candidate_count': len(candidates), 'paired_cells': len(cells), 'terminal_slots': slots,
        'comparison': parent['comparison'], 'new_molecular_calls': 0})
    return candidates, cells, old_ids, helper, base, world


def import_solvent(ev, candidate, metal, cell):
    sid, x = candidate['source_id'], candidate['_x']
    raw_pin = read(cell['point'])['ALPB'] if 'point' in cell else cell['ALPB']
    raw = read(raw_pin)
    for p in raw['artifacts'].values():
        verify(p)
    start, execution = read(raw['artifacts']['START.json']), read(raw['artifacts']['EXECUTION.json'])
    require(execution['returncode'] == 0 and execution['command'] == start['command']
            and same_bytes(pin(start['command'][0]), ev.cfg['ALPB_binary'])
            and start['command'][1:] == [str(verify(start['xyz'])), str(verify(start['charges']))],
            'Module2 ALPB native binary/command failed or differs')
    symbols, rx = ev.admit.xyz(start['xyz'])
    q = read(ev.tasks[(sid, 'baseline', metal)]['global_solvent_charges'])['charges_e']
    require(q == read(ev.tasks[(sid, 'expanded', metal)]['global_solvent_charges'])['charges_e']
            and [float(v) for v in verify(start['charges']).read_text().split()] == q
            and symbols == [a['element'] for a in ev.partitions[(sid, 'baseline')][1]['physical_atoms'][:-1]] + [metal]
            and ev.np.array_equal(rx, x) and raw['rounding_bounds_hartree'] == [0., 0.]
            and abs(raw['charge_e'] + 2) <= 1e-5 and raw['total_hartree'] == cell['native_energies']['ALPB_Eh']
            and math.fsum(raw['components_hartree'].values()) == raw['total_hartree'],
            'Module2 ALPB exact physical state or native energy differs')
    grad = ev.np.loadtxt(verify(raw['artifacts']['gradient.dat']))
    require(grad.shape == (256, 3) and ev.np.isfinite(grad).all(), 'Module2 native ALPB gradient unavailable')
    value = {'energy_hartree': raw['total_hartree'], 'gradient_hartree_per_bohr': grad.tolist(),
             'components_hartree': raw['components_hartree'], 'record': raw_pin, 'source': 'exact_admitted_Module2_cell'}
    key = ev.solvent_key(sid, metal, x, q)
    if key in ev.solvent_cache:
        require(all(ev.solvent_cache[key][k] == value[k] for k in
                    ('energy_hartree', 'gradient_hartree_per_bohr', 'components_hartree')),
                'Conflicting exact global solvent evidence')
    else:
        ev.solvent_cache[key] = value
    return {'candidate_id': candidate['candidate_id'], 'source_id': sid, 'metal': metal,
            'ALPB': raw_pin, 'solvent_cache_key': key, 'status': 'complete'}


def cell_from_point(slot, p, reused):
    slot.update(status='complete', point=p['provenance']['point'], cache_key=p['provenance']['cache_key'], reused=reused,
                native_energies=p['native_energies'], component_energies_kcal_mol=p['component_energies_kcal_mol'],
                total_energy_kcal_mol=p['total_energy_kcal_mol'])


def menu(base, sid, region, ids, cells, origin_id):
    missing = [{'candidate_id': c, 'metal': m} for c in ids for m in METALS
               if cells.get((c, region, m), {}).get('status') != 'complete']
    value = {'source_id': sid, 'region': region, 'candidate_ids': ids, 'required_cells': 2 * len(ids),
             'status': 'unavailable', 'missing_cells': missing}
    if missing:
        return value
    works = {c: {m: base.difference(cells[(c, region, m)], cells[(origin_id, region, m)]) for m in METALS} for c in ids}
    selected = {m: min(ids, key=lambda c: works[c][m]['composite_kcal_mol']) for m in METALS}
    selected_work = {m: works[selected[m]][m] for m in METALS}
    value.update(status='available', selected=selected,
                 exact_ties={m: [c for c in ids if works[c][m]['composite_kcal_mol'] == selected_work[m]['composite_kcal_mol']] for m in METALS},
                 selected_work_from_original_kcal_mol=selected_work,
                 common_accommodation_kcal_mol={k: math.fsum(selected_work[m][k] for m in METALS) / 2 for k in selected_work['La']},
                 differential_accommodation_Dy_minus_La_kcal_mol={k: selected_work['Dy'][k] - selected_work['La'][k] for k in selected_work['La']})
    return value


def prepare(cfg, out, result):
    require(cfg['protocol'] == PROTOCOL and same_bytes(cfg['implementation'], pin(__file__)), 'Pool source/protocol differs')
    verify(cfg['specification'])
    verify(cfg['plan'])
    base, ev = qualified_hybrid(cfg, out, result)
    for region in REGIONS:
        a, b = [ev.partitions[(s, region)][1] for s in SOURCES]
        identity = lambda p: [(x['id'], x['element'], x['kind'], x.get('resname')) for x in p['physical_atoms']]
        field = lambda p: [(x['id'], x['source_index'], x['element'], x['charge_e']) for x in p['field']]
        require(identity(a) == identity(b) and a['qm_source_indices'] == b['qm_source_indices']
                and a['caps'] == b['caps'] and a['qm_charge'] == b['qm_charge'] and field(a) == field(b),
                'Same-region source chemical/partition inventory differs')
        for metal in METALS:
            require(read(ev.tasks[('7CCO', region, metal)]['global_solvent_charges'])['charges_e'] ==
                    read(ev.tasks[('7CCN', region, metal)]['global_solvent_charges'])['charges_e'],
                    'Same-region source fixed charge state differs')
    candidates, old_cells, old_ids, whole_helper, old_base, world = module2_pool(cfg, out, ev, result)
    origin_ids = {}
    for c in candidates:
        sid, x = c['source_id'], c['_x']
        physical = read(c['physical_source'])
        require(len(physical) == 256 and ev.np.array_equal([a['xyz_A'] for a in physical], x), 'Pool physical coordinates differ')
        for region in REGIONS:
            atoms = ev.partitions[(sid, region)][1]['physical_atoms']
            require([{k: v for k, v in a.items() if k != 'xyz_A'} for a in physical] ==
                    [{k: v for k, v in a.items() if k != 'xyz_A'} for a in atoms], 'Hybrid/Module2 physical inventory differs')
        if ev.np.array_equal(x, ev.origins[(sid, 'baseline')]):
            require(sid not in origin_ids, 'Duplicated source origin')
            origin_ids[sid] = c['candidate_id']
    require(set(origin_ids) == set(SOURCES), 'Both original source geometries required')
    slots, cells = [], {}
    for sid in SOURCES:
        for region in REGIONS:
            for metal in METALS:
                for c in candidates:
                    if c['source_id'] != sid:
                        continue
                    slot = {'candidate_id': c['candidate_id'], 'source_id': sid, 'region': region, 'metal': metal,
                            'status': 'unavailable_not_attempted', 'physical_source': c['physical_source']}
                    key = ev.key(sid, region, metal, c['_x'])
                    if key in ev.cache:
                        require(ev.cache[key] is not None, 'Prior exact hybrid failure cannot be retried')
                        cell_from_point(slot, ev.cache[key], True)
                    slots.append(slot)
                    cells[(c['candidate_id'], region, metal)] = slot
    require(all(cells[(origin_ids[s], r, m)]['status'] == 'complete' for s in SOURCES for r in REGIONS for m in METALS),
            'Eight exact origins were not reused')
    missing = sum(r['status'] != 'complete' for r in slots)
    require(len(slots) == 4 * len(candidates) and missing <= 4 * len(candidates) - 8 <= 96, 'Hybrid pool denominator exceeds declaration')
    test_cells = dict(cells)
    removed = (origin_ids['7CCO'], 'baseline', 'La')
    test_cells.pop(removed)
    test = menu(base, '7CCO', 'baseline', [origin_ids['7CCO']], test_cells, origin_ids['7CCO'])
    require(test['status'] == 'unavailable' and test['missing_cells'] == [{'candidate_id': origin_ids['7CCO'], 'metal': 'La'}],
            'Actual required-cell omission test failed')
    require(all(ev.cost[k] == 0 for k in ('composites_attempted', 'MACE_calls_attempted', 'GFN_calls_attempted',
                'MM_configurations_attempted', 'MM_group_queries_attempted', 'ALPB_calls_attempted')) and ev.native is None,
            'Preflight made a molecular call')
    result.update(candidates=[{k: v for k, v in c.items() if k != '_x'} for c in candidates],
                  inherited_candidate_ids=old_ids, origin_candidate_ids=origin_ids, declared_geometries=len(candidates),
                  declared_hybrid_cells=len(slots), cells=slots, compatible_hybrid_cells_reused=len(slots) - missing,
                  maximum_new_composites=missing, identity_sha256=ev.identity_hash,
                  component_call_ceiling={'MACE': missing, 'GFN2': 2 * missing, 'MM_configurations': missing,
                                          'MM_group_queries': 4 * missing, 'ALPB': 0})
    solvent_replays = [{'candidate_id': c['candidate_id'], 'source_id': c['source_id'], 'metal': m,
                        'status': 'unavailable_not_admitted'} for c in candidates for m in METALS]
    result['required_prior_solvent_cells'] = solvent_replays
    by_id = {c['candidate_id']: c for c in candidates}
    for admission in solvent_replays:
        cid, metal = admission['candidate_id'], admission['metal']
        try:
            admission.update(import_solvent(ev, by_id[cid], metal, old_cells[(cid, metal)]))
        except Exception as exc:
            admission.update(status='unavailable_prior_evidence', failure=repr(exc))
            raise
    result['preflight'] = write(out / 'PREFLIGHT.json', {'status': 'passed', 'new_molecular_calls': 0,
        'model_initialized': False, 'hybrid_parent_replay': result['hybrid_parent_replay'],
        'same_region_cross_source_chemical_identity': True,
        'module2_replay': result['module2_replay'], 'candidates': result['candidates'],
        'declared_cells': slots, 'maximum_new_composites': missing, 'component_call_ceiling': result['component_call_ceiling'],
        'solvent_replays': solvent_replays, 'missing_actual_cell_test': {'test_only': True, 'removed': list(removed), 'result': test}})
    return base, ev, candidates, cells, old_cells, whole_helper


def summarize(base, candidates, cells, old_cells, whole_helper, result):
    origins = result['origin_candidate_ids']
    menus, inherited, geometry = [], [], []
    for region in REGIONS:
        for sid in SOURCES:
            ids = [c['candidate_id'] for c in candidates if c['source_id'] == sid]
            menus.append(menu(base, sid, region, ids, cells, origins[sid]))
            inherited.append(menu(base, sid, region, [c for c in ids if c in result['inherited_candidate_ids']], cells, origins[sid]))
    for c in candidates:
        cid, sid = c['candidate_id'], c['source_id']
        row = {'candidate_id': cid, 'source_id': sid, 'status': 'unavailable', 'physical_source': c['physical_source']}
        geometry.append(row)
        required = [cells[(cid, r, m)] for r in REGIONS for m in METALS]
        if any(p['status'] != 'complete' for p in required):
            continue
        work = {r: {m: base.difference(cells[(cid, r, m)], cells[(origins[sid], r, m)]) for m in METALS} for r in REGIONS}
        whole = {m: whole_helper.difference(old_cells[(cid, m)]['native_energies'], old_cells[(origins[sid], m)]['native_energies']) for m in METALS}
        shift = {m: {k: work['expanded'][m][k] - work['baseline'][m][k] for k in work['baseline'][m]} for m in METALS}
        row.update(status='available', within_region_work_from_original_kcal_mol=work,
                   whole_model_work_from_original_components_kcal_mol=whole,
                   same_geometry_expanded_minus_baseline_work_kcal_mol=shift,
                   same_geometry_hybrid_minus_whole_total_work_kcal_mol={r: {m: work[r][m]['composite_kcal_mol'] - math.fsum(whole[m].values()) for m in METALS} for r in REGIONS},
                   region_differential_accommodation_Dy_minus_La_kcal_mol={r: {k: work[r]['Dy'][k] - work[r]['La'][k] for k in work[r]['La']} for r in REGIONS},
                   same_geometry_partition_effect_on_differential_work_kcal_mol={k: shift['Dy'][k] - shift['La'][k] for k in shift['La']})
    spreads = []
    for region in REGIONS:
        own = {m['source_id']: m for m in menus if m['region'] == region}
        row = {'region': region, 'status': 'unavailable', 'definition': '(Dy7CCN-Dy7CCO) - (La7CCN-La7CCO), within one region'}
        spreads.append(row)
        if not all(m['status'] == 'available' for m in own.values()):
            continue
        values = {}
        for scope in ('original', 'selected'):
            ids = {s: {m: origins[s] if scope == 'original' else own[s]['selected'][m] for m in METALS} for s in SOURCES}
            strains = {m: base.difference(cells[(ids['7CCN'][m], region, m)], cells[(ids['7CCO'][m], region, m)]) for m in METALS}
            differential = {k: strains['Dy'][k] - strains['La'][k] for k in strains['La']}
            values[scope] = {'same_metal_7CCN_minus_7CCO_kcal_mol': strains, 'differential_source_spread_kcal_mol': differential,
                             'absolute_total_differential_source_spread_kcal_mol': abs(differential['composite_kcal_mol'])}
        row.update(status='available', **values)
    parent_selected = []
    geometry_by_id = {r['candidate_id']: r for r in geometry}
    for old_menu in result['module2_comparison']['source_menus']:
        for metal in METALS:
            cid = old_menu['selected'][metal]
            g = geometry_by_id[cid]
            row = {'source_id': old_menu['source_scope'], 'metal': metal, 'candidate_id': cid, 'status': g['status']}
            if g['status'] == 'available':
                row['hybrid_work_at_same_Module2_selected_geometry_kcal_mol'] = {r: g['within_region_work_from_original_kcal_mol'][r][metal] for r in REGIONS}
            parent_selected.append(row)
    result['comparison'] = {'complete_source_menus': menus, 'inherited_22_candidate_submenus': inherited,
                            'source_spreads': spreads, 'same_geometry_work_comparisons': geometry,
                            'Module2_selected_geometries': parent_selected,
                            'no_cross_source_pool': True, 'no_unlike_region_absolute_energy_comparison': True}
    return all(m['status'] == 'available' for m in menus)


def report(out, result):
    lines = ['# Complete LBT pool under both hybrid partitions', '', 'Status: ' + result['status'] + '.', '',
             '| Source | Region | Availability | La work | Dy work | Dy−La accommodation |',
             '| --- | --- | --- | ---: | ---: | ---: |']
    for m in result.get('comparison', {}).get('complete_source_menus', []):
        if m['status'] == 'available':
            w = m['selected_work_from_original_kcal_mol']
            lines.append('| %s | %s | available | %.6f | %.6f | %.6f |' % (m['source_id'], m['region'],
                w['La']['composite_kcal_mol'], w['Dy']['composite_kcal_mol'], m['differential_accommodation_Dy_minus_La_kcal_mol']['composite_kcal_mol']))
        else:
            lines.append('| %s | %s | unavailable | — | — | — |' % (m['source_id'], m['region']))
    lines += ['', 'Works are kcal/mol relative to the original geometry under the same source, metal and region.',
              'All admitted Module2 candidates are retained. Per-geometry partition effects compare these works,',
              'never absolute energies between partitions. Source spreads and component terms remain in RESULT.json.',
              'No optimization, affinity, population, electronic-state accuracy or convergence claim is made.', '']
    (out / 'REPORT.md').write_text('\n'.join(lines))
    result['report'] = pin(out / 'REPORT.md')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage', choices=('preflight', 'run'), required=True)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick = time.monotonic()
    initial_usage = resource.getrusage(resource.RUSAGE_SELF)
    result = {'protocol': PROTOCOL, 'status': 'unavailable_prior_evidence', 'stage': a.stage,
              'implementation': pin(__file__), 'config': pin(a.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
              'new_DFT_calls': 0, 'optimizer_steps': 0, 'optimization_qualified': False}
    try:
        cfg = json.loads(a.config.read_text())
        base, ev, candidates, cells, old_cells, whole_helper = prepare(cfg, out, result)
        if a.stage == 'preflight':
            result['status'] = 'complete'
            result['preflight_only'] = True
        else:
            result['status'] = 'unavailable_backend'
            by_id = {c['candidate_id']: c for c in candidates}
            for slot in result['cells']:
                if slot['status'] == 'complete':
                    continue
                c = by_id[slot['candidate_id']]
                before, request_count = copy.deepcopy(ev.cost), len(ev.requests)
                slot['status'] = 'unavailable_attempted'
                try:
                    value = ev.evaluate(slot['source_id'], slot['region'], slot['metal'], c['_x'])
                    cell_from_point(slot, value, False)
                except Exception as exc:
                    slot.update(failure=repr(exc), traceback=traceback.format_exc())
                    break
                finally:
                    slot['cost_delta'] = {k: ev.cost[k] - before.get(k, 0) for k in ev.cost}
                    slot['attempted_point'] = ev.requests[-1] if len(ev.requests) > request_count else None
            limit = result['maximum_new_composites']
            require(ev.cost['composites_attempted'] <= limit and ev.cost['MACE_calls_attempted'] <= limit
                    and ev.cost['GFN_calls_attempted'] <= 2 * limit and ev.cost['MM_configurations_attempted'] <= limit
                    and ev.cost['MM_group_queries_attempted'] <= 4 * limit and ev.cost['ALPB_calls_attempted'] == 0,
                    'Declared call ceiling or exact solvent reuse violated')
            if summarize(base, candidates, cells, old_cells, whole_helper, result):
                result['status'] = 'complete'
        verify(result['config'])
        verify(cfg['specification'])
    except Exception as exc:
        if result['status'] == 'complete':
            result['status'] = 'execution_failure'
        result.update(failure=repr(exc), traceback=traceback.format_exc())
    ev = result.pop('_evaluator', None)
    if ev is not None:
        result.update(actual_cost=ev.cost, requests=ev.requests, points=ev.points)
        if ev.native is not None:
            result.update(native_MACE_calls=ev.native.calls, native_MACE_device=ev.native.device,
                          native_MACE_versions=ev.native.versions,
                          peak_CUDA_allocated_bytes=ev.native.torch.cuda.max_memory_allocated())
    result['elapsed_seconds'] = time.monotonic() - tick
    usage = resource.getrusage(resource.RUSAGE_SELF)
    result.update(peak_host_RSS_KiB=usage.ru_maxrss,
                  observed_parent_CPU_seconds=usage.ru_utime + usage.ru_stime - initial_usage.ru_utime - initial_usage.ru_stime,
                  application_allocated_CPU_seconds=result['elapsed_seconds'] * int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
    report(out, result)
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
