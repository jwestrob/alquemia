#!/usr/bin/env python3
"""Saved six-search mechanics comparison. No molecular backend or submission."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import time
import traceback

PROTOCOL = 'lanm-six-search-broad-response-saved-comparison-v1'
STATES = (('LBT', '7CCO', 'La'), ('LBT', '7CCO', 'Dy'),
          ('LBT', '7CCN', 'La'), ('LBT', '7CCN', 'Dy'),
          ('Mex', '8FNS', 'La'), ('Mex', '8FNS', 'Dy'))
COMPONENTS = ('MACE_vacuum', 'GFN2_field_difference', 'MM_ledger', 'fixed_FF_charge_ALPB')
EV = 23.06054783061903
HA = 627.509474
ARITHMETIC_TOLERANCE = 1e-7
RESOLUTION = 1.0


def require(test, message):
    if not test:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(record):
    require(isinstance(record, dict) and set(record) >= {'path', 'sha256'}, 'Expected exact file pin')
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: '+record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def same_pin(a, b):
    verify(a); verify(b)
    return a['sha256'] == b['sha256']


def close(a, b, label, tolerance=ARITHMETIC_TOLERANCE):
    require(math.isfinite(a) and math.isfinite(b) and abs(a-b) <= tolerance,
            label+': '+repr((a, b)))


def parts(value):
    require(set(value) >= set(COMPONENTS), 'Missing energy component')
    result = {k: float(value[k]) for k in COMPONENTS}
    require(all(math.isfinite(v) for v in result.values()), 'Nonfinite component')
    return result


def native_work(new, old, system):
    """Same stable difference equations as existing MEX_RESPONSE.work and LBT helper."""
    result = {'MACE_vacuum': (new['MACE_eV']-old['MACE_eV'])*EV,
              'fixed_FF_charge_ALPB': (new['ALPB_Eh']-old['ALPB_Eh'])*HA}
    if system == 'Mex':
        require(set(new['MM_groups_kcal_mol']) == set(old['MM_groups_kcal_mol']), 'MM groups changed')
        result['GFN2_field_difference'] = math.fsum((new['GFN2_embedded_Eh']-old['GFN2_embedded_Eh'],
            -(new['GFN2_isolated_Eh']-old['GFN2_isolated_Eh'])))*HA
        result['MM_ledger'] = math.fsum(new['MM_groups_kcal_mol'][k]-old['MM_groups_kcal_mol'][k]
                                       for k in old['MM_groups_kcal_mol'])
    else:
        require(new['field_Eh'] == old['field_Eh'] == 0.0 and
                new['MM_kcal_mol'] == old['MM_kcal_mol'] == 0.0, 'Whole-LBT recipe changed')
        result.update(GFN2_field_difference=0.0, MM_ledger=0.0)
    return parts(result)


def point(provenance, system):
    record = provenance if system == 'LBT' else provenance['point']
    value = read(record)
    require(value['status'] == 'complete', 'Saved point is unavailable')
    if system == 'Mex':
        request = read(provenance['cartesian_request'])
        require(same_pin(request['point'], record) and request['cache_key'] == provenance['cache_key'] and
                request['identity_sha256'] == provenance['identity_sha256'] and request['state'] == provenance['state']
                and request['native_energies'] == value['native_energies'] and
                same_pin(request['gradient'], value['assembled_gradient']), 'Mex point/request identity differs')
        gradient = read(value['assembled_gradient'])['gradient_kcal_mol_A']
    else:
        import numpy as np
        gradient = np.load(verify(value['gradient_kcal_mol_A']), allow_pickle=False).tolist()
    return value, gradient


def force_metrics(gradient, atom_count):
    require(len(gradient) == atom_count and all(len(row) == 3 for row in gradient), 'Full physical gradient size differs')
    squared = [math.fsum(float(v)*float(v) for v in row) for row in gradient]
    require(all(math.isfinite(v) for v in squared), 'Nonfinite physical gradient')
    maximum = math.sqrt(max(squared)); rms = math.sqrt(math.fsum(squared)/atom_count)
    return {'maximum_atom_kcal_mol_A': maximum, 'RMS_atom_vector_kcal_mol_A': rms,
            'L2_kcal_mol_A': math.sqrt(math.fsum(squared)), 'converged': maximum <= .10 and rms <= .03}


def inspect_search(slot, baseline, scout):
    system, sid, metal = slot['system'], slot['source_id'], slot['metal']
    row = dict(system=system, source_id=sid, metal=metal, status='unavailable', result=slot['result'],
               initial_force=None, final_force=None, component_work_kcal_mol=None, total_work_kcal_mol=None,
               converged=False, physical_state_admitted=False, prior_attempts=slot['attempt_results'])
    if slot['result'] is None:
        row['reason'] = 'No result supplied'; return row
    parent = read(slot['result']); cfg = read(parent['config']); verify(parent['implementation'])
    require(same_pin(parent['implementation'], cfg['implementation']), 'Run implementation/config differs')
    row.update(execution_status=parent['status'], failure=parent.get('failure'), slurm_job_id=parent.get('slurm_job_id'))
    if system == 'LBT':
        require(same_pin(parent['parent_pool_result'], baseline['old_lbt_pool']), 'LBT run old pool differs')
        rows = parent.get('rows', [])
        if len(rows) != 1:
            row['reason'] = 'Run lacks its single terminal search'; return row
        require((rows[0]['source_id'], rows[0]['metal']) == (sid, metal), 'LBT terminal state differs')
        terminal = rows[0]; search_pin = terminal['search']
    else:
        require(parent['metal'] == metal, 'Mex run target differs')
        ecfg = read(cfg['evaluator_config'])
        require(same_pin(ecfg['pool_result'], baseline['old_mex_pool']), 'Mex run old pool differs')
        if 'search' not in parent:
            row['reason'] = 'Run lacks its terminal search'; return row
        terminal = parent; search_pin = parent['search']
    search = read(search_pin)
    row.update(search=search_pin, termination=search['termination'], accepted_steps=search['accepted_steps'],
               distinct_evaluations=search['distinct_evaluations'], limits=search['limits'],
               physical_state_admitted=bool(search['physical_state_admitted']),
               search_status=search['status'], failure=search.get('failure'))
    require(search['accepted_steps'] <= 80 and search['distinct_evaluations'] <= 200 and
            search['settings']['maximum_force'] == .10 and
            search['settings']['RMS_atom_force'] == .03, 'Search budget or full-force thresholds changed')
    initial, gi = point(search['start_point'], system); final, gf = point(search['terminal_point'], system)
    require(initial['metal'] == final['metal'] == metal, 'Point target differs')
    if system == 'LBT':
        require(initial['source_id'] == final['source_id'] == sid and initial['recipe'] == final['recipe'],
                'LBT point recipe/source changed')
        old = baseline['old_rows'][sid]['new_menu']
        require(initial['native_energies'] == old['selected_native_energies'][metal], 'LBT start is not old selected minimum')
        if sid == '7CCO' and metal == 'La':
            require(scout is not None, 'Original scout required for cumulative accounting')
            ss = read(scout['rows'][0]['search'])
            require(ss['start_point'] == search['start_point'], 'Scout original start differs')
    else:
        old = baseline['mex']['common5']['selections'][metal]['point']
        require(initial['native_energies'] == read(old)['native_energies'], 'Mex start is not old selected minimum')
        require(search['start_point']['identity_sha256'] == search['terminal_point']['identity_sha256'] and
                search['start_point']['state'] == search['terminal_point']['state'], 'Mex state/recipe changed')
    atom_count = 256 if system == 'LBT' else 2049
    fi, ff = force_metrics(gi, atom_count), force_metrics(gf, atom_count)
    for k in ('maximum_atom_kcal_mol_A', 'RMS_atom_vector_kcal_mol_A', 'L2_kcal_mol_A'):
        close(ff[k], search['final_full_force'][k], 'Final full-force replay '+k, 1e-10+1e-12*abs(ff[k]))
    require(ff['converged'] == search['converged'], 'Full-force convergence differs')
    work = native_work(final['native_energies'], initial['native_energies'], system)
    total = math.fsum(work.values()); close(total, search['work_from_start']['composite_kcal_mol'], 'Own search work replay')
    reported = search['work_from_start']
    if system == 'LBT':
        close(work['MACE_vacuum'], reported['MACE_kcal_mol'], 'LBT MACE work')
        close(work['fixed_FF_charge_ALPB'], reported['fixed_ALPB_kcal_mol'], 'LBT ALPB work')
    else:
        for k in COMPONENTS:
            close(work[k], reported[k], 'Mex component work '+k)
    admitted = parent['status'] == search['status'] == 'complete' and search['physical_state_admitted']
    if system == 'LBT':
        admitted = admitted and terminal['status'] == 'complete' and terminal['physical_state_admitted']
    row.update(status='admitted_finite_response' if admitted else 'unavailable',
               initial_force=fi, final_force=ff, converged=ff['converged'],
               component_work_kcal_mol=work, total_work_kcal_mol=total,
               initial_point=search['start_point'], terminal_point=search['terminal_point'],
               final_geometry=search['final_geometry'])
    if not admitted:
        row['reason'] = 'Execution or bound-state admission failed; terminal work retained descriptively'
    return row


def pool_rows(cfg, baseline, searches):
    result = []
    lbt = None if cfg['lbt_pool_result'] is None else read(cfg['lbt_pool_result'])
    mex = None if cfg['mex_pool_result'] is None else read(cfg['mex_pool_result'])
    for system, sid in (('LBT', '7CCO'), ('LBT', '7CCN'), ('Mex', '8FNS')):
        parent = lbt if system == 'LBT' else mex
        row = {'system': system, 'source_id': sid, 'status': 'unavailable', 'selected_work_kcal_mol': None,
               'selected_work_components_kcal_mol': None, 'Dy_minus_La_change_kcal_mol': None}
        result.append(row)
        if parent is None:
            row['reason'] = 'No paired-pool result supplied'; continue
        if parent.get('stage') != 'pool':
            row['reason'] = 'Result is not the completed paired-pool stage'; continue
        pconfig = read(parent['config']); verify(parent['implementation'])
        require(same_pin(parent['implementation'], pconfig['implementation']), 'Pool implementation differs')
        required = [s for s in searches if s['system'] == system and s['source_id'] == sid]
        supplied_pins = pconfig['run_results'] if system == 'LBT' else pconfig['search_results']
        configured = [s['result'] for s in cfg['searches'] if s['system'] == system]
        require(all(a is not None for a in configured) and len(configured) == len(supplied_pins)
                and all(same_pin(a, b) for a, b in zip(configured, supplied_pins)), 'Pool/search result lineage differs')
        if system == 'LBT':
            require(same_pin(parent['old_parent_pool'], cfg['old_lbt_pool']), 'LBT pool baseline differs')
            old = next(x for x in parent['old_source_baselines'] if x['source_scope'] == sid)
            expected = baseline['old_rows'][sid]['new_menu']
            for key in ('selected', 'selected_native_energies', 'LBT_Dy_minus_La_components_kcal_mol'):
                require(old[key] == expected[key], 'Old LBT menu replay differs: '+key)
            menus = parent.get('comparison', {}).get('source_menus', [])
            menu = next((x for x in menus if x['source_scope'] == sid), None)
        else:
            require(same_pin(parent['old_pool_result'], cfg['old_mex_pool']), 'Mex pool baseline differs')
            old = parent['old_common5']; expected = baseline['mex']['common5']
            require(all(old[k] == v for k, v in expected.items()), 'Old Mex common5 replay differs')
            menu = parent.get('common7')
        row['menu'] = menu; row['old_menu'] = old
        if parent.get('status') not in ('complete', 'unavailable'):
            row['reason'] = 'Pool execution failed; retained menu is descriptive only'; continue
        if menu is None or menu.get('status') != 'available' or any(s['status'] != 'admitted_finite_response' for s in required):
            row['reason'] = 'A declared search or required paired state is unavailable'; continue
        components = {}; totals = {}
        for metal in ('La', 'Dy'):
            if system == 'LBT':
                work = native_work(menu['selected_native_energies'][metal], old['selected_native_energies'][metal], system)
                reported = menu['selected_work_from_old_minimum_components_kcal_mol'][metal]
                expected_total = menu['selected_work_from_old_minimum_kcal_mol'][metal]
            else:
                new_native = read(menu['selections'][metal]['point'])['native_energies']
                old_native = read(old['selections'][metal]['point'])['native_energies']
                work = native_work(new_native, old_native, system)
                reported = parent['selected_work_from_common5'][metal]
                expected_total = reported['composite_kcal_mol']
            for key in COMPONENTS:
                close(work[key], reported[key], 'Selected pool work '+system+'/'+metal+'/'+key)
            components[metal] = work; totals[metal] = math.fsum(work.values())
            close(totals[metal], expected_total, 'Selected pool work total')
            require(totals[metal] <= ARITHMETIC_TOLERANCE, 'Expanded pool excludes the old minimum')
        change = math.fsum((totals['Dy'], -totals['La']))
        expected_change = (menu['contrast_change_from_old_kcal_mol'] if system == 'LBT'
                           else parent['common5_to_common7_Dy_minus_La_change_kcal_mol'])
        close(change, expected_change, 'Paired response change')
        row.update(status='available', selected_work_components_kcal_mol=components,
                   selected_work_kcal_mol=totals, Dy_minus_La_change_kcal_mol=change,
                   all_searches_converged=all(s['converged'] for s in required))
    return result


def conditional_rows(baseline, pools):
    mex = next(p for p in pools if p['system'] == 'Mex'); rows = []
    for sid in ('7CCO', '7CCN'):
        old = baseline['old_rows'][sid]; lbt = next(p for p in pools if p['system'] == 'LBT' and p['source_id'] == sid)
        row = {'source_scope': sid, 'status': 'unavailable', 'old_R_kcal_mol': old['new_R_kcal_mol'],
               'old_R_components_kcal_mol': old['new_R_components_kcal_mol'],
               'new_R_kcal_mol': None, 'new_R_components_kcal_mol': None,
               'unmatched_scaffold_Hamiltonian': True, 'qualified_preference_prediction': False}
        rows.append(row)
        close(math.fsum(parts(row['old_R_components_kcal_mol']).values()), row['old_R_kcal_mol'], 'Old R closure')
        if mex['status'] != 'available' or lbt['status'] != 'available':
            row['reason'] = 'Missing admitted source-local paired menu'; continue
        changes = {k: math.fsum((mex['selected_work_components_kcal_mol']['Dy'][k],
            -mex['selected_work_components_kcal_mol']['La'][k], -lbt['selected_work_components_kcal_mol']['Dy'][k],
            lbt['selected_work_components_kcal_mol']['La'][k])) for k in COMPONENTS}
        new_parts = {k: math.fsum((row['old_R_components_kcal_mol'][k], changes[k])) for k in COMPONENTS}
        change = math.fsum(changes.values()); new = math.fsum(new_parts.values())
        close(new, math.fsum((row['old_R_kcal_mol'], change)), 'New conditional R closure')
        row.update(status='available', new_R_kcal_mol=new, new_R_components_kcal_mol=new_parts,
            R_change_components_kcal_mol=changes, R_change_kcal_mol=change,
            absolute_change_exceeds_1_kcal_mol=abs(change)>RESOLUTION,
            within_1_kcal_mol_resolution=abs(new)<=RESOLUTION,
            all_searches_converged=mex['all_searches_converged'] and lbt['all_searches_converged'])
    return rows


def cost_ledger(cfg):
    artifacts = []; seen = set()
    def add(label, record):
        if record is None or record['sha256'] in seen:
            return
        seen.add(record['sha256']); value = read(record)
        artifacts.append({'label': label, 'result': record, 'status': value.get('status'),
            'slurm_job_id': value.get('slurm_job_id'), 'elapsed_seconds': value.get('elapsed_seconds'),
            'process_CPU_seconds': value.get('process_CPU_seconds'),
            'actual_costs': value.get('actual_costs'), 'actual_cost': value.get('actual_cost'),
            'counters': value.get('counters'), 'allocation': value.get('allocation'),
            'application_allocated_CPU_seconds': value.get('application_allocated_CPU_seconds'),
            'application_allocated_GPU_seconds': value.get('application_allocated_GPU_seconds'),
            'peak_host_RSS_KiB': value.get('peak_host_RSS_KiB'),
            'failure': value.get('failure')})
    add('original LBT7CCO La scout', cfg['scout_result'])
    for slot in cfg['searches']:
        label = '/'.join(slot[k] for k in ('system', 'source_id', 'metal'))
        add(label+' final bound result', slot['result'])
        for i, record in enumerate(slot['attempt_results']):
            add(label+' prior attempt '+str(i+1), record)
    add('LBT paired pool', cfg['lbt_pool_result']); add('Mex paired pool', cfg['mex_pool_result'])
    jobs = {}; receipt_sources = []
    for record in cfg['cost_receipts']:
        receipt = read(record); receipt_sources.append(record)
        require('records' in receipt, 'Expected immutable completion-watcher receipt snapshot')
        for jid, values in receipt['records'].items():
            if jid in jobs:
                require(jobs[jid]['recorded'] == values, 'Conflicting scheduler accounting for '+jid)
                continue
            elapsed = int(values['elapsed_seconds']); cpus = int(values['allocated_cpus'])
            jobs[jid] = {'recorded': values, 'scheduler_elapsed_seconds': elapsed, 'allocated_CPUs': cpus,
                         'allocated_CPU_seconds': elapsed*cpus, 'receipt': record}
    owned = {str(r['slurm_job_id']) for r in artifacts if r['slurm_job_id'] is not None}
    require(len(owned) == sum(r['slurm_job_id'] is not None for r in artifacts),
            'Multiple distinct top-level result receipts for the same allocation; bind one authoritative result')
    require(set(jobs) <= owned, 'Cost receipt includes an undeclared job; do not count watcher/unrelated jobs')
    calls = {}
    mappings = {'composites': ('new_composite_evaluations', 'new_composite_evaluations_attempted'),
                'MACE': ('MACE_calls_attempted', 'new_MACE_calls_attempted'),
                'GFN2': ('new_GFN2_calls_attempted',), 'MM_configurations': ('new_MM_configurations_attempted',),
                'MM_group_queries': ('new_MM_group_queries_attempted',), 'ALPB': ('ALPB_calls_attempted', 'new_ALPB_calls_attempted')}
    for name, aliases in mappings.items():
        reports = []
        for artifact in artifacts:
            counter = artifact['actual_costs'] or artifact['actual_cost'] or artifact['counters']
            entries = counter if isinstance(counter, list) else [counter]
            values = []
            for entry in entries:
                require(entry is None or isinstance(entry, dict), 'Unknown native counter schema')
                matches = [] if entry is None else [entry[k] for k in aliases if k in entry]
                require(len(matches) <= 1, 'Ambiguous per-allocation counter aliases')
                values.append(None if not matches else matches[0])
            reports.append(sum(values) if values and all(v is not None for v in values) else None)
        calls[name] = {'reported_subtotal': sum(v for v in reports if v is not None),
                       'artifacts_with_counter': sum(v is not None for v in reports),
                       'artifacts_without_counter': sum(v is None for v in reports)}
    return {'artifacts': artifacts, 'scheduler_jobs': jobs, 'scheduler_receipts': receipt_sources,
            'jobs_without_scheduler_receipt': sorted(owned-set(jobs)),
            'allocated_CPU_seconds_known_jobs': sum(x['allocated_CPU_seconds'] for x in jobs.values()),
            'component_call_reported_subtotals': calls, 'utilization_measured': False,
            'notes': ['Original scout and each exact result are counted once.',
                'Cumulative optimizer counters are not summed with per-allocation component counters.',
                'Absent counters remain unknown; reported subtotals are not complete-cost claims.',
                'Observed process CPU excludes unreported children; scheduler allocation is not utilization.']}


def table(path, rows, fields):
    with Path(path).open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(row.get(k), separators=(',', ':')) if isinstance(row.get(k), (dict, list))
                             else row.get(k) for k in fields})


def run(cfg, out, result):
    require(same_pin(cfg['implementation'], pin(__file__)), 'Comparison implementation changed')
    verify(cfg['plan']); verify(cfg['specification'])
    require([(s['system'], s['source_id'], s['metal']) for s in cfg['searches']] == list(STATES), 'Initial six denominator changed')
    old = read(cfg['old_conditional_result']); lbt = read(cfg['old_lbt_pool']); mex = read(cfg['old_mex_pool'])
    oldcfg = read(old['config'])
    require(old['status'] == lbt['status'] == mex['status'] == 'complete' and
            old['slurm_job_id'] == '12769' and lbt['slurm_job_id'] == '12767' and mex['slurm_job_id'] == '12731',
            'Exact completed restricted baseline required')
    require(same_pin(oldcfg['response_pool'], cfg['old_lbt_pool']), 'Conditional baseline LBT lineage differs')
    oldrows = {x['source_scope']: x for x in old['rows'] if x['source_scope'] in ('7CCO', '7CCN')}
    require(set(oldrows) == {'7CCO', '7CCN'}, 'Old source-local denominator differs')
    for row in oldrows.values():
        require(row['status'] == 'available' and row['Mex_selected_unchanged'] == mex['common5']['selections'], 'Old conditional Mex selections differ')
    baseline = dict(old_rows=oldrows, mex=mex, old_lbt_pool=cfg['old_lbt_pool'], old_mex_pool=cfg['old_mex_pool'])
    scout = None if cfg['scout_result'] is None else read(cfg['scout_result'])
    searches = result['searches']
    for i, slot in enumerate(cfg['searches']):
        searches[i] = inspect_search(slot, baseline, scout)
    pools = pool_rows(cfg, baseline, searches); result['source_local_pools'] = pools
    conditional = conditional_rows(baseline, pools)
    old_signed = oldrows['7CCN']['new_R_kcal_mol']-oldrows['7CCO']['new_R_kcal_mol']
    spread = {'status': 'unavailable', 'old_R_7CCN_minus_7CCO_kcal_mol': old_signed, 'old_absolute_kcal_mol': abs(old_signed)}
    if all(r['status'] == 'available' for r in conditional):
        new_signed = conditional[1]['new_R_kcal_mol']-conditional[0]['new_R_kcal_mol']
        spread.update(status='available', new_R_7CCN_minus_7CCO_kcal_mol=new_signed,
                      new_absolute_kcal_mol=abs(new_signed), absolute_spread_change_kcal_mol=abs(new_signed)-abs(old_signed))
    masked = copy.deepcopy(pools)
    next(p for p in masked if p['system'] == 'Mex')['status'] = 'unavailable'
    test = conditional_rows(baseline, masked)
    require(all(r['status'] == 'unavailable' and r['new_R_kcal_mol'] is None for r in test), 'Unavailable-pool suppression failed')
    electronic = read(cfg['old_electronic_substitution'])
    require(electronic['status'] == 'complete' and electronic['slurm_job_id'] == '12823', 'Original electronic diagnostic unavailable')
    costs = cost_ledger(cfg)
    result.update(status='complete', declared_search_slots=6, searches=searches, source_local_pools=pools,
        conditional_exchange=conditional, source_spread=spread,
        admitted_finite_searches=sum(s['status'] == 'admitted_finite_response' for s in searches),
        full_force_converged_searches=sum(s['status'] == 'admitted_finite_response' and s['converged'] for s in searches),
        available_conditional_rows=sum(r['status'] == 'available' for r in conditional),
        tests={'six_slots_preserved': len(searches) == 6, 'old_selected_components_replayed': True,
               'native_work_and_force_replayed_for_available_searches': True,
               'test_only_masked_Mex_menu_suppresses_both_new_R_rows': True},
        separate_electronic_sensitivity_reference={'result': cfg['old_electronic_substitution'],
            'source_local_spreads': electronic.get('source_local_spreads'),
            'used_in_new_response_arithmetic': False,
            'scope': 'Earlier LBT pool and changed electronic model. Model sensitivity, not GFN accuracy.'},
        interpretation={'unmatched_scaffold_Hamiltonian': True, 'binding_free_energy': False,
            'qualified_preference_prediction': False, 'experimental_validation': False,
            'LBT_7CCN_conditioning': 'Lu', 'Mex_state': 'EF3 La/Dy target; three Nd spectators and retained waters',
            'motion_test': 'Broad versus restricted under each unchanged per-system energy model',
            'resolution_kcal_mol': RESOLUTION, 'resolution_is_uncertainty': False,
            'source_cross_pool_selection': False})
    write(out/'COSTS.json', costs); result['costs'] = pin(out/'COSTS.json')
    table(out/'SEARCHES.tsv', searches, ['system', 'source_id', 'metal', 'status', 'termination', 'accepted_steps',
        'distinct_evaluations', 'total_work_kcal_mol', 'component_work_kcal_mol', 'initial_force', 'final_force', 'converged', 'reason'])
    table(out/'POOLS.tsv', pools, ['system', 'source_id', 'status', 'selected_work_kcal_mol',
        'selected_work_components_kcal_mol', 'Dy_minus_La_change_kcal_mol', 'all_searches_converged', 'reason'])
    table(out/'CONDITIONAL.tsv', conditional, ['source_scope', 'status', 'old_R_kcal_mol', 'new_R_kcal_mol',
        'R_change_kcal_mol', 'R_change_components_kcal_mol', 'all_searches_converged', 'reason'])
    lines = ['# Broad-coordinate response under unchanged per-system energies', '',
        f"{result['admitted_finite_searches']}/6 finite searches admitted; {result['full_force_converged_searches']}/6 satisfy full-force convergence.", '',
        '| System/source | Metal | Additional work (kcal/mol) | Final max/RMS force | Outcome |',
        '| --- | --- | ---: | --- | --- |']
    for s in searches:
        f = s['final_force']; force = 'unavailable' if f is None else f"{f['maximum_atom_kcal_mol_A']:.3f}/{f['RMS_atom_vector_kcal_mol_A']:.3f}"
        work = 'unavailable' if s['total_work_kcal_mol'] is None else f"{s['total_work_kcal_mol']:.6f}"
        lines.append(f"| {s['system']} {s['source_id']} | {s['metal']} | {work} | {force} | {s.get('termination', s['status'])} |")
    lines += ['', 'Forces are kcal/mol/Å. Finite energy lowering is not a stationary minimum or accuracy test.', '',
        '## Conditional exchange retained for continuity', '',
        '**The scaffold Hamiltonians remain unmatched.** LBT uses whole-peptide MACE; Mex uses local electronics plus MM. These are not qualified preference predictions.', '',
        '| LBT source | Old R | New R | Change (kcal/mol) |', '| --- | ---: | ---: | ---: |']
    for r in conditional:
        new = 'unavailable' if r['new_R_kcal_mol'] is None else f"{r['new_R_kcal_mol']:.6f}"
        change = 'unavailable' if r.get('R_change_kcal_mol') is None else f"{r['R_change_kcal_mol']:.6f}"
        lines.append(f"| {r['source_scope']} | {r['old_R_kcal_mol']:.6f} | {new} | {change} |")
    lines += ['', '7CCN is Lu-conditioned. Mex retains three Nd spectators. No cross-source minimum, affinity estimate or experimental validation is implied.', '',
        '## Decision', '', 'Use the measured component work and residual forces to diagnose missing mechanics. Module 3 must test matched local-electronic/MM representation; independent references must test accuracy. Energy lowering, a smaller source gap or a desired sign alone is insufficient.', '',
        f"Known scheduler allocation cost: {costs['allocated_CPU_seconds_known_jobs']} CPU-s; missing receipts: {', '.join(costs['jobs_without_scheduler_receipt']) or 'none'}. Component calls, measured wall time and missing counters are preserved separately in COSTS.json.", '']
    (out/'REPORT.md').write_text('\n'.join(lines))
    result['report'] = pin(out/'REPORT.md')


def bind(record):
    if record is None:
        return None
    if isinstance(record, str):
        return pin(record)
    verify(record); return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage', choices=('template', 'bind', 'compare'), required=True)
    p.add_argument('--repo', type=Path); p.add_argument('--config', type=Path); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.stage == 'template':
        require(args.repo is not None, '--repo required for template')
        repo = args.repo.resolve(); base = repo/'workspaces/lanm_coupled_response_20261016'
        config = {'protocol_id': PROTOCOL, 'plan': pin(repo/'diagnostics/lanm_coupled_mechanics_20261017/PLAN.md'),
            'old_conditional_result': pin(base/'monomer_matched_response_comparison_v1/execution/RESULT.json'),
            'old_lbt_pool': pin(base/'lbt_matched_response_pool_v1/execution/RESULT.json'),
            'old_mex_pool': pin(base/'mex_continuation_pool_v1/execution/RESULT.json'),
            'old_electronic_substitution': pin(base/'electronic_substitution_v2/execution/RESULT.json'),
            'searches': [dict(system=s, source_id=sid, metal=m, result=None, attempt_results=[]) for s, sid, m in STATES],
            'scout_result': None, 'lbt_pool_result': None, 'mex_pool_result': None, 'cost_receipts': []}
        write(args.output, config); return 0
    require(args.config is not None, '--config required')
    cfg = json.loads(args.config.read_text())
    if args.stage == 'bind':
        cfg['implementation'] = pin(__file__); cfg['specification'] = pin(Path(__file__).with_name('SPEC.md'))
        for key in ('plan', 'old_conditional_result', 'old_lbt_pool', 'old_mex_pool', 'old_electronic_substitution',
                    'scout_result', 'lbt_pool_result', 'mex_pool_result'):
            cfg[key] = bind(cfg[key])
        for slot in cfg['searches']:
            slot['result'] = bind(slot['result']); slot['attempt_results'] = [bind(r) for r in slot['attempt_results']]
        cfg['cost_receipts'] = [bind(r) for r in cfg['cost_receipts']]
        write(args.output, cfg); return 0
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Saved numerical analysis requires Slurm')
    require(cfg['protocol_id'] == PROTOCOL, 'Protocol differs')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic(); cpu = time.process_time()
    result = {'protocol_id': PROTOCOL, 'status': 'execution_failure', 'implementation': pin(__file__),
              'config': pin(args.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'new_molecular_calls': 0,
              'new_geometries': 0, 'optimizer_steps': 0, 'declared_search_slots': 6,
              'searches': [dict(system=s, source_id=sid, metal=m, status='unavailable_not_reviewed') for s, sid, m in STATES]}
    try:
        run(cfg, out, result); verify(result['config'])
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick, process_CPU_seconds=time.process_time()-cpu)
    write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
