#!/usr/bin/env python3
"""Saved whole/hybrid/response accounting; no molecular evaluation path."""
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

PROTOCOL = 'lanm_hybrid_model_and_response_saved_comparison_v1'
STATES = (('7CCO', 'La'), ('7CCO', 'Dy'), ('7CCN', 'La'), ('7CCN', 'Dy'))
SOURCES = ('7CCO', '7CCN')
METALS = ('La', 'Dy')
STAGES = ('Module2_whole', 'hybrid_same_Module2_selection', 'hybrid_parent_pool', 'hybrid_response_pool')
FROZEN = {'module2_implementation': 'f9cdb6a4f6e36bde97f5250d01f24cf6bb51af9fca37f8e600fce77b3198939a',
          'response_implementation': '05b1c412feebfc72048405c8723d426168872325ce32d956732f5c898c028981'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        return {'path': str(p), 'sha256': hashlib.file_digest(f, 'sha256').hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Changed pin: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def same(a, b):
    verify(a)
    verify(b)
    return a['sha256'] == b['sha256']


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


def replay_module2(cfg, out):
    parent = read(cfg['module2_result'])
    pcfg = read(parent['config'])
    require(parent['status'] == 'complete' and parent['slurm_job_id'] == '12982'
            and same(parent['implementation'], cfg['module2_implementation'])
            and same(pcfg['implementation'], cfg['module2_implementation']), 'Complete fixed Module 2 result required')
    legacy = module('hybrid_comparison_module2', pcfg['implementation'])
    replay = {'searches': [None] * 6}
    d = out / 'module2_replay'
    d.mkdir()
    legacy.run(pcfg, d, replay)
    for key in ('searches', 'source_local_pools', 'conditional_exchange', 'source_spread',
                'admitted_finite_searches', 'full_force_converged_searches', 'available_conditional_rows'):
        require(replay[key] == parent[key], 'Module 2 saved comparison differs: ' + key)
    require(parent['admitted_finite_searches'] == 6 and parent['available_conditional_rows'] == 2,
            'Complete six-state Module 2 comparison required')
    return legacy, parent, pcfg


def replay_mex(legacy_cfg, out):
    """Use saved admission/summary only; never call the pool's scoring path."""
    parent = read(legacy_cfg['mex_pool_result'])
    cfg = read(parent['config'])
    require(parent['status'] == 'complete' and parent['stage'] == 'pool'
            and parent['declared_geometries'] == 7 and parent['declared_cells'] == 14
            and same(parent['implementation'], cfg['implementation']), 'Fixed complete Mex common7 required')
    pool = module('hybrid_comparison_mex_pool', cfg['implementation'])
    ecfg = read(cfg['evaluator_config'])
    require(same(ecfg['implementation'], cfg['evaluator_implementation']), 'Mex evaluator lineage differs')
    em = module('hybrid_comparison_mex_evaluator', ecfg['implementation'])
    sm = module('hybrid_comparison_mex_search', cfg['search_implementation'])
    core = module('hybrid_comparison_mex_core', cfg['optimizer'])
    d = out / 'mex_replay'
    d.mkdir()
    ev = em.MexEvaluator(ecfg, d / 'saved_inputs')
    guard, _ = sm.build_guard(ev, core, d)
    oldpool = module('hybrid_comparison_mex_old_pool', ev.cfg['pool_implementation'])
    matrix, by_coordinate, old, _ = pool.baseline(ev, oldpool, guard, d)
    slots, failures = pool.admit_searches(cfg, ev, sm, core, guard, matrix, by_coordinate)
    require(all(s['status'] == 'admitted' for s in slots) and not failures, 'Mex terminal/failure admission differs')
    for request in parent['requests']:
        raw = read(request)
        pool.import_request(ev, raw['state']['metal'], request)
    require(len(matrix) == len(parent['matrix']) == 7, 'Mex common7 denominator differs')
    for row, saved in zip(matrix, parent['matrix']):
        require(row['coordinate_sha256'] == saved['coordinate_sha256'] and row['aliases'] == saved['aliases'],
                'Mex exact coordinate/alias menu differs')
        for metal in METALS:
            cell = pool.cached_cell(ev, metal, row['_x'])
            require(cell['native_energies'] == saved['cells'][metal]['native_energies']
                    and same(cell['point'], saved['cells'][metal]['point']), 'Mex raw paired state differs')
            row['cells'][metal] = cell
        pool.update_work(ev, row)
    summary = pool.summarize(oldpool, ev, matrix, pool.ALIASES, slots)
    require(summary == parent['common7'] and old == parent['old_common5'] and ev.new_evaluations == 0,
            'Mex common7 selections differ or replay made a call')
    selected = {m: read(summary['selections'][m]['point'])['native_energies'] for m in METALS}
    return parent, selected, ev.identity, ev.state_metadata, ev.old.origins['build']


def replay_response(cfg, out, result, legacy_cfg):
    parent = read(cfg['response_pool_result'])
    rcfg = read(parent['config'])
    require(parent['status'] == 'complete' and parent['stage'] == 'pool' and parent['new_pool_available']
            and parent['region'] == rcfg['region'] == 'baseline'
            and same(parent['implementation'], cfg['response_implementation'])
            and same(rcfg['implementation'], cfg['response_implementation'])
            and same(rcfg['plan'], cfg['plan'])
            and same(rcfg['parent_pool_result'], cfg['hybrid_parent_result']), 'Complete baseline response pool required')
    response = module('hybrid_comparison_response', rcfg['implementation'])
    d = out / 'response_replay'
    d.mkdir()
    scratch = {}
    loaded = response.prepare(rcfg, d, scratch)
    ev, base, pool = loaded['ev'], loaded['base'], loaded['pool']
    hp = loaded['parent']
    hcfg = read(hp['config'])
    require(same(hcfg['module2_pool_result'], legacy_cfg['lbt_pool_result']), 'Hybrid pool uses a different Module 2 pool')
    result['hybrid_parent_comparison'] = hp['comparison']
    require(len(rcfg['run_results']) == len(parent['terminal_slots']) == 4, 'Four response search slots required')
    terminals = []
    for state, record, slot in zip(STATES, rcfg['run_results'], parent['terminal_slots']):
        require((slot['source_id'], slot['metal']) == state and slot['status'] == 'admitted'
                and same(slot['result'], record), 'Response terminal slot differs')
        _, search, _, point = response.saved_search(record, rcfg, loaded, state, ('run',))
        require(point['provenance'] == slot['point'] and same(slot['search'], read(record)['search'])
                and search['final_full_force'] == slot['full_force']
                and search['converged'] == slot['converged'] and search['termination'] == slot['termination'],
                'Response terminal force or provenance differs')
        terminals.append(dict(slot, initial_full_force=next(s['initial_full_force'] for s in scratch['starts']
                              if (s['source_id'], s['metal']) == state), work_from_start=search['work_from_start']))
    candidates = parent['candidates']
    require(26 <= len(candidates) <= 30 and parent['declared_geometries'] == len(candidates)
            and parent['declared_cells'] == 2 * len(candidates), 'Full response pool denominator differs')
    inherited = [{k: v for k, v in c.items() if k != '_x'} for c in loaded['candidates']]
    for new, old in zip(candidates[:26], inherited):
        require(all(new[k] == v for k, v in old.items() if k != 'aliases')
                and new['aliases'][:len(old['aliases'])] == old['aliases'], 'Inherited candidate changed or missing')
    coordinates, seen = {}, set()
    for index, c in enumerate(candidates):
        x = ev.np.asarray([a['xyz_A'] for a in read(c['physical_source'])], dtype='<f8')
        key = (c['source_id'], x.tobytes())
        require(c['source_id'] in SOURCES and key not in seen and c['candidate_id'] not in coordinates
                and x.shape == (256, 3), 'Response source coordinate identity/dedup differs')
        if index < 26:
            # Historical candidates have two pinned hash schemas. Their original
            # physical pin and full state have already replayed in prepare().
            require(ev.np.array_equal(x, loaded['candidates'][index]['_x']), 'Inherited exact coordinates differ')
        else:
            require(hashlib.sha256(x.tobytes()).hexdigest() == c['physical_coordinate_bytes_sha256'],
                    'New terminal coordinate bytes differ')
        seen.add(key)
        coordinates[c['candidate_id']] = (c['source_id'], x)
    for slot in terminals:
        sid, x = coordinates[slot['candidate_id']]
        point = loaded['qualifier'].restore(ev, slot['point']['point'])
        candidate = next(c for c in candidates if c['candidate_id'] == slot['candidate_id'])
        require(sid == slot['source_id'] and ev.np.array_equal(x, point['positions_A'])
                and sid + ':hybrid_baseline_terminal_' + slot['metal'] in candidate['aliases'],
                'Terminal alias or physical geometry differs')
    require({c['candidate_id'] for c in candidates[26:]} <= {s['candidate_id'] for s in terminals},
            'New candidate is not an admitted response terminal')
    expected = [(c['candidate_id'], 'baseline', m) for c in candidates for m in METALS]
    actual = [(r['candidate_id'], r['region'], r['metal']) for r in parent['paired_cells']]
    require(actual == expected, 'Complete ordered response paired denominator differs')
    cells = {}
    for row, key in zip(parent['paired_cells'], actual):
        require(row['status'] == 'complete', 'Required response paired cell unavailable')
        point = loaded['qualifier'].restore(ev, row['point'])
        sid, x = coordinates[key[0]]
        require((point['source_id'], point['region'], point['metal']) == (sid, key[1], key[2])
                and ev.np.array_equal(point['positions_A'], x)
                and all(row[k] == point[k] for k in ('native_energies', 'component_energies_kcal_mol', 'total_energy_kcal_mol'))
                and row['cache_key'] == point['provenance']['cache_key'], 'Response native point binding differs')
        cells[key] = row
    menus = {}
    for sid, saved in zip(SOURCES, parent['comparison']['source_menus']):
        require(saved['source_id'] == sid and saved['region'] == 'baseline', 'Response source menu order differs')
        ids = [c['candidate_id'] for c in candidates if c['source_id'] == sid]
        menu = pool.menu(base, sid, 'baseline', ids, cells, hp['origin_candidate_ids'][sid])
        old = next(m for m in hp['comparison']['complete_source_menus'] if (m['source_id'], m['region']) == (sid, 'baseline'))
        require(menu == saved['new'] and old == saved['old'] and menu['status'] == 'available', 'Response selection replay differs')
        works = {m: base.difference(cells[(menu['selected'][m], 'baseline', m)], loaded['cells'][(old['selected'][m], 'baseline', m)]) for m in METALS}
        require(works == saved['additional_selected_work_from_parent_minimum_kcal_mol']
                and {k: works['Dy'][k] - works['La'][k] for k in works['La']} == saved['change_in_Dy_minus_La_kcal_mol'],
                'Response selected component work differs')
        menus[sid] = menu
    masked = dict(cells)
    first = expected[0]
    del masked[first]
    sid = candidates[0]['source_id']
    test = pool.menu(base, sid, 'baseline', menus[sid]['candidate_ids'], masked, hp['origin_candidate_ids'][sid])
    require(test['status'] == 'unavailable' and len(test['missing_cells']) == 1, 'Real missing-cell suppression failed')
    require(ev.cost['composites_attempted'] == 0 and ev.native is None, 'Saved response replay made a molecular call')
    result.update(response_terminal_slots=terminals, response_pool_comparison=parent['comparison'],
                  missing_cell_test={'removed_actual_cell': list(first), 'test_only': True, 'result': test})
    return loaded, parent, cells, menus


def compare(cfg, out, result):
    require(cfg['protocol'] == PROTOCOL and same(cfg['implementation'], pin(__file__)), 'Comparator source/protocol differs')
    verify(cfg['specification'])
    verify(cfg['plan'])
    for key, sha in FROZEN.items():
        require(verify(cfg[key]) and cfg[key]['sha256'] == sha, 'Frozen comparison dependency differs: ' + key)
    missing = [k for k in ('module2_result', 'hybrid_parent_result', 'response_pool_result') if cfg[k] is None]
    if missing:
        result.update(status='unbound', missing_required_parents=missing)
        return
    legacy, old, oldcfg = replay_module2(cfg, out)
    mex, mex_native, mex_identity, mex_states, mex_build = replay_mex(oldcfg, out)
    loaded, response, final_cells, final_menus = replay_response(cfg, out, result, oldcfg)
    base, hp, ev = loaded['base'], loaded['parent'], loaded['ev']
    require(base.EV == legacy.EV and base.HA == legacy.HA, 'Native energy conversion differs')
    require(mex_identity['MACE_config']['model'] == ev.mcfg['model']
            and mex_identity['GFN_runtime'] == ev.runtime, 'Mex/LBT local electronic model or standalone runtime differs')
    # The published native functional is shared, but scope assertions produce
    # deliberately different binaries for 2049 and 256 physical atoms.
    require(same(mex_identity['solvent_binary'], mex_build['binary'])
            and same(ev.alpb['build']['binary'], ev.cfg['ALPB_binary']), 'Native solvent binary provenance differs')
    mex_fortran = verify(mex_build['fortran']).read_text()
    marker = "if (n /= 2049) error stop 'declared Mex physical scope requires 2049 atoms'"
    require(mex_fortran.count(marker) == 1 and verify(ev.alpb['build']['fortran']).read_text()
            == mex_fortran.replace(marker, "if (n /= 256) error stop 'declared LBT physical scope requires 256 atoms'")
            and mex_identity['solvent_config']['build_inputs'] == ev.alpb['build']['build_inputs'],
            'Native solvent functional/build differs beyond the declared atom-count assertion')
    hybrid_states = []
    for sid, metal in STATES:
        state = ev.states[(sid, 'baseline', metal)]
        mstate = mex_states[metal]
        require(state['metal'] == mstate['metal'] == metal and state['oxidation_state'] == mstate['oxidation_state'] == 3
                and state['MACE_multiplicity'] == mstate['physical_multiplicity'] == (1 if metal == 'La' else 6)
                and state['GFN_uhf'] == mstate['GFN2_state']['uhf'] == 0
                and state['GFN_effective_multiplicity'] == mstate['GFN2_state']['multiplicity'] == 1,
                'One-target physical/effective electronic exchange state differs')
        hybrid_states.append(state)
    mex_delta = legacy.native_work(mex_native['Dy'], mex_native['La'], 'Mex')
    rows = []
    for sid in SOURCES:
        oldrow = next(r for r in old['conditional_exchange'] if r['source_scope'] == sid)
        whole = next(r for r in hp['module2_comparison']['source_menus'] if r['source_scope'] == sid)
        hybrid = next(r for r in hp['comparison']['complete_source_menus'] if (r['source_id'], r['region']) == (sid, 'baseline'))
        choices = {'Module2_whole': whole['selected'], 'hybrid_same_Module2_selection': whole['selected'],
                   'hybrid_parent_pool': hybrid['selected'], 'hybrid_response_pool': final_menus[sid]['selected']}
        stage_rows = []
        for stage in STAGES:
            if stage == 'Module2_whole':
                delta = legacy.native_work(whole['selected_native_energies']['Dy'], whole['selected_native_energies']['La'], 'LBT')
            else:
                cells = final_cells if stage == 'hybrid_response_pool' else loaded['cells']
                delta = base.difference(cells[(choices[stage]['Dy'], 'baseline', 'Dy')], cells[(choices[stage]['La'], 'baseline', 'La')])
            parts = {k: math.fsum((mex_delta[k], -delta[k])) for k in legacy.COMPONENTS}
            total = math.fsum(parts.values())
            if stage == 'Module2_whole':
                legacy.close(total, oldrow['new_R_kcal_mol'], 'Unchanged Module 2 R')
                for k in legacy.COMPONENTS:
                    legacy.close(parts[k], oldrow['new_R_components_kcal_mol'][k], 'Unchanged Module 2 R component')
            stage_rows.append({'stage': stage, 'selected': choices[stage], 'R_components_kcal_mol': parts,
                               'R_kcal_mol': total, 'abs_R_le_1_kcal_mol': abs(total) <= legacy.RESOLUTION})
        changes = []
        labels = ('same_geometry_Hamiltonian_change', 'existing_pool_reselection', 'additional_hybrid_response')
        for label, a, b in zip(labels, stage_rows, stage_rows[1:]):
            parts = {k: math.fsum((b['R_components_kcal_mol'][k], -a['R_components_kcal_mol'][k])) for k in legacy.COMPONENTS}
            changes.append({'term': label, 'R_change_components_kcal_mol': parts, 'R_change_kcal_mol': math.fsum(parts.values())})
        for k in legacy.COMPONENTS:
            legacy.close(math.fsum(c['R_change_components_kcal_mol'][k] for c in changes),
                         stage_rows[-1]['R_components_kcal_mol'][k] - stage_rows[0]['R_components_kcal_mol'][k], 'Ordered change closure')
        rows.append({'source_id': sid, 'status': 'available', 'stages': stage_rows, 'ordered_changes': changes})
    spreads = []
    for i, stage in enumerate(STAGES):
        parts = {k: rows[1]['stages'][i]['R_components_kcal_mol'][k] - rows[0]['stages'][i]['R_components_kcal_mol'][k] for k in legacy.COMPONENTS}
        signed = math.fsum(parts.values())
        spreads.append({'stage': stage, 'R_7CCN_minus_7CCO_components_kcal_mol': parts,
                        'R_7CCN_minus_7CCO_kcal_mol': signed, 'absolute_kcal_mol': abs(signed)})
    scouts, previous_runs = {}, {}
    for slot in response['terminal_slots']:
        run_cfg = read(read(slot['result'])['config'])
        for key, ledger in (('scout_result', scouts), ('resume_result', previous_runs)):
            record = run_cfg.get(key)
            if record is not None:
                ledger[record['sha256']] = {'result': record, 'actual_cost': read(record)['actual_cost']}
    result.update(status='complete', rows=rows, source_spreads=spreads, Module2_unchanged=old['conditional_exchange'],
                  Module2_searches=old['searches'], Mex_fixed_common7=mex['common7'],
                  Mex_terminal_slots=mex['terminal_slots'], Mex_pool_result=oldcfg['mex_pool_result'],
                  electronic_identity={'model': ev.mcfg['model'], 'GFN_runtime': ev.runtime,
                    'LBT_ALPB_binary': ev.cfg['ALPB_binary'], 'Mex_ALPB_binary': mex_build['binary'],
                    'ALPB_functional_differs_only_in_atom_count_assertion': True,
                    'Mex_states': mex_states, 'LBT_states': hybrid_states},
                  parent_costs={'Module2_cost_ledger': old['costs'], 'hybrid_parent': hp['actual_cost'],
                                'response_pool': response['actual_cost'],
                                'response_scouts_separate': list(scouts.values()),
                                'response_prior_complete_runs_separate': list(previous_runs.values()),
                                'allocation_costs_not_summed_with_cumulative_search_counts': True,
                                'response_searches': [{'result': s['result'], 'actual_cost': read(s['result'])['actual_cost']} for s in response['terminal_slots']]},
                  interpretation={'binding_free_energy': False, 'qualified_preference_prediction': False,
                    'cross_source_pool': False, 'cross_region_absolute_energy_subtraction': False,
                    'LBT_7CCN_conditioning': 'Lu', 'MACE_physical_multiplicities': {'La': 1, 'Dy': 6},
                    'GFN2_effective_multiplicity': 1, 'frozen_f_approximation_retained': True,
                    'resolution_kcal_mol': legacy.RESOLUTION, 'resolution_is_uncertainty': False})


def report(out, result):
    lines = ['# Saved scaffold-model and response comparison', '', 'Status: ' + result['status'] + '.', '']
    if result['status'] == 'complete':
        lines += ['R = (Mex Dy−La) − (LBT Dy−La), kcal/mol. Mex common7 stays fixed.', '',
                  '| LBT source | Module 2 whole | Hybrid, same coordinates | Hybrid parent pool | Hybrid response pool |',
                  '| --- | ---: | ---: | ---: | ---: |']
        for row in result['rows']:
            lines.append('| ' + row['source_id'] + ' | ' + ' | '.join('%.6f' % r['R_kcal_mol'] for r in row['stages']) + ' |')
        lines += ['', '| Source | Same-coordinate model change | Existing-pool reselection | Added hybrid response |',
                  '| --- | ---: | ---: | ---: |']
        for row in result['rows']:
            lines.append('| ' + row['source_id'] + ' | ' + ' | '.join('%.6f' % r['R_change_kcal_mol'] for r in row['ordered_changes']) + ' |')
        lines += ['', '| Stage | R(7CCN) − R(7CCO) | Absolute source difference |', '| --- | ---: | ---: |']
        for row in result['source_spreads']:
            lines.append('| %s | %.6f | %.6f |' % (row['stage'], row['R_7CCN_minus_7CCO_kcal_mol'], row['absolute_kcal_mol']))
        lines += ['', '| Source/metal | Final max/RMS force | Converged | Termination |', '| --- | ---: | --- | --- |']
        for s in result['response_terminal_slots']:
            f = s['full_force']
            lines.append('| %s/%s | %.4f / %.4f | %s | %s |' % (s['source_id'], s['metal'],
                f['maximum_atom_kcal_mol_A'], f['RMS_atom_vector_kcal_mol_A'], s['converged'], s['termination']))
        lines += ['', 'Forces are kcal/mol/Å. Ordered component changes and both source differences remain in RESULT.json.',
                  'Finite lowering is not stationary convergence or electronic accuracy. The Lu-conditioned source,',
                  'different inventories/assays and frozen-f GFN2 approximation prevent an affinity claim.',
                  'Expanded-region sensitivity remains separate in the full parent comparison; it does not select the primary model.']
    else:
        lines += ['A complete final comparison is unavailable. No partial parent or earlier model substitutes for missing evidence.']
    (out / 'REPORT.md').write_text('\n'.join(lines) + '\n')
    result['report'] = pin(out / 'REPORT.md')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Saved numerical comparison requires Slurm')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick, cpu = time.monotonic(), time.process_time()
    result = {'protocol': PROTOCOL, 'status': 'execution_failure', 'implementation': pin(__file__),
              'config': pin(a.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
              'new_molecular_calls': 0, 'new_geometries': 0, 'optimizer_steps': 0,
              'declared_sources': list(SOURCES), 'declared_stages': list(STAGES), 'declared_response_states': [list(s) for s in STATES]}
    try:
        cfg = json.loads(a.config.read_text())
        compare(cfg, out, result)
        verify(result['config'])
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic() - tick, process_CPU_seconds=time.process_time() - cpu)
    report(out, result)
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
