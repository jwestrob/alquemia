"""Frozen-band source-conditioning robustness summaries; no fitting or execution."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import json
from pathlib import Path
import statistics

from affordable_common import HA_TO_KCAL, InvalidArtifact, cache_key, energy, read_json, record, verify, write_new, xyz
from compact_solvation_compare import COMPOSITE_PROTOCOL, mix_pair, native_endpoint, pinned_manifest

SCHEMA = 'PQQ_fold_conditioning_comparison_v1'
REPS = ('core', 'context')
COMPONENTS = ('native_R_model_kcal_mol', 'GFN2_vacuum_R_kcal_mol',
              'GFN2_ALPB_R_kcal_mol', 'solvation_delta_R_kcal_mol', 'composite_R_model_kcal_mol')


@lru_cache(maxsize=32)
def original_tasks(path, sha256):
    from run_orca_task_manifest import load_manifest_tasks
    pin = {'path': path, 'sha256': sha256}
    m, tasks = load_manifest_tasks(verify(pin))
    return m, {t['task_id']: t for t in tasks}


def original_completed(pin, task_id):
    """Same receipt check as native collector, with immutable task-map caching."""
    from run_orca_task_manifest import _completed_attempt_is_valid
    m, tasks = original_tasks(pin['path'], pin['sha256']); task = tasks[task_id]
    out = task['output']; receipt = Path(str(out) + '.execution.json')
    if (not out.exists() or not receipt.exists() or not _completed_attempt_is_valid(
            receipt, out, manifest_sha256=pin['sha256'], task=task,
            runner_identity=m['execution_policy']['task_runner'],
            runtime_renderer_identity=m['execution_policy']['runtime_renderer'])):
        return None
    return {'energy_hartree': energy(out), 'output': record(out), 'receipt': record(receipt),
            'manifest': pin, 'task_id': task_id}


def checked_scoped_transfer(row, source, inventory_pin, current_manifest):
    """Validate original receipts without relabeling explicitly reused endpoints."""
    from compact_solvation import METHOD, input_text
    m = pinned_manifest(current_manifest['path'], current_manifest['sha256'])
    if (m['inventory'] != inventory_pin or m['method_id'] != METHOD
            or m['protocol_id'] != COMPOSITE_PROTOCOL or row['status'] not in ('complete', 'unavailable')):
        raise InvalidArtifact('current transfer scope/method differs')
    complete_media = [medium for medium in ('vacuum', 'alpb') if row['endpoints'][medium]['status'] == 'complete']
    if (len(complete_media) == 2) != (row['status'] == 'complete'):
        raise InvalidArtifact('transfer availability disagrees with its endpoints')
    reused = {}
    for medium in complete_media:
        matches = [t for t in m['all_tasks'] if (t['case_id'], t['representation'], t['metal'], t['medium'], t['solver']) ==
                   (row['case_id'], row['representation'], row['metal'], medium, row['solver'])]
        if len(matches) != 1: raise InvalidArtifact('missing/ambiguous declared transfer task')
        task = matches[0]; e = row['endpoints'][medium]
        if task['source_endpoint'] != source or e['status'] != 'complete':
            raise InvalidArtifact('transfer source or status differs')
        body = input_text(source['charge'], source['multiplicity'], medium, row['solver'])
        scientific_key = cache_key({'xyz_sha256': source['xyz']['sha256'], 'charge': source['charge'],
            'multiplicity': source['multiplicity'], 'medium': medium, 'solver': row['solver'],
            'method': METHOD, 'input_body': body, 'orca': m['orca']})
        if (task['scientific_key'] != scientific_key or verify(task['input']).read_text() != body
                or verify(task['xyz']).read_bytes() != verify(source['xyz']).read_bytes()):
            raise InvalidArtifact('current calculation does not match its scientific key')
        actual = original_completed(e['manifest'], e['task_id'])
        if actual is None or any(actual[k] != e[k] for k in ('energy_hartree', 'output', 'receipt', 'manifest', 'task_id')):
            raise InvalidArtifact('original transfer receipt does not replay')
        if task['task_id'] in m['reused']:
            if actual != m['reused'][task['task_id']]:
                raise InvalidArtifact('reused receipt is not the declared original endpoint')
            original = pinned_manifest(e['manifest']['path'], e['manifest']['sha256'])
            oldtask = next(t for t in original['all_tasks'] if t['task_id'] == e['task_id'])
            if (original['method_id'] != METHOD or original['protocol_id'] != COMPOSITE_PROTOCOL
                    or oldtask['scientific_key'] != scientific_key
                    or (oldtask['charge'], oldtask['multiplicity'], oldtask['metal'], oldtask['medium'], oldtask['solver']) !=
                       (source['charge'], source['multiplicity'], row['metal'], medium, row['solver'])
                    or verify(oldtask['xyz']).read_bytes() != verify(source['xyz']).read_bytes()):
                raise InvalidArtifact('reused original physical state/scientific key differs')
            reused[medium] = {'current_task_id': task['task_id'], 'original_manifest': e['manifest'],
                              'original_task_id': e['task_id'], 'scientific_key': scientific_key}
        elif e['manifest'] != current_manifest or e['task_id'] != task['task_id']:
            raise InvalidArtifact('receipt belongs to an undeclared reuse')
        if row[medium + '_hartree'] != actual['energy_hartree']:
            raise InvalidArtifact('transfer energy does not match original receipt')
    expected_transfer = row['alpb_hartree'] - row['vacuum_hartree'] if len(complete_media) == 2 else None
    if row['delta_solv_hartree'] != expected_transfer:
        raise InvalidArtifact('transfer sign/algebra differs')
    return reused


def source_groups(source):
    cases = source['cases']; ids = [c['case_id'] for c in cases]
    if len(cases) != 250 or len(set(ids)) != 250:
        raise InvalidArtifact('frozen scope requires 250 distinct sources')
    grouped = defaultdict(list)
    for case in cases:
        grouped[case['root_case_id']].append(case)
    if len(grouped) != 25:
        raise InvalidArtifact('frozen scope requires 25 parent proteins')
    for root, members in grouped.items():
        ca = [m for m in members if m['source_conditioning_metal'] == 'Ca']
        la = [m for m in members if m['source_conditioning_metal'] == 'La']
        canonical = [m for m in members if m['canonical_coordinate_match']]
        if (len(members), len(ca), len(la), len(canonical)) != (10, 5, 5, 1):
            raise InvalidArtifact('incorrect source arms/replay count: ' + root)
        if canonical[0]['source_conditioning_metal'] != 'La':
            raise InvalidArtifact('canonical replay must belong to the La source arm')
        if any(m['primary_evaluation_pool'] == m['canonical_coordinate_match'] for m in members):
            raise InvalidArtifact('primary pool must exclude exactly the canonical source')
        if len({(m['expected_class'], m['biological_group']) for m in members}) != 1:
            raise InvalidArtifact('within-protein source labels/groups differ')
    return dict(grouped)


def bands(reference):
    if reference['protocol_id'] != COMPOSITE_PROTOCOL or reference['solver'] != 'native':
        raise InvalidArtifact('reference must be the existing native composite calibration')
    result = {}
    for rep in REPS:
        composite = reference['calibration'][rep]['bands']
        if not composite or composite['representation'] != rep:
            raise InvalidArtifact('missing compatible frozen bands')
        native = reference['native_PQQ_bands'][rep]
        result[rep] = {
            'composite': {'Ca_max': composite['Ca_supported_max_R_model_kcal_mol'],
                          'La_min': composite['La_supported_min_R_model_kcal_mol']},
            'native': {'Ca_max': native['U_max_Ca_kcal_mol'] if rep == 'core' else native['Ca_supported_max_R_kcal_mol'],
                       'La_min': native['L_min_La_kcal_mol'] if rep == 'core' else native['La_supported_min_R_kcal_mol']}}
        if any(b['Ca_max'] >= b['La_min'] for b in result[rep].values()):
            raise InvalidArtifact('frozen class bands overlap')
    return result


def decision(value, frozen):
    if value is None: return 'unavailable'
    if value <= frozen['Ca_max']: return 'Ca-supported'
    if value >= frozen['La_min']: return 'La-supported'
    return 'inconclusive'


def pool_summary(members, index, frozen):
    """Every declared member is required independently for each component."""
    rows = [index[m['case_id']] for m in members]
    result = {'declared_members': [m['case_id'] for m in members], 'denominator': len(members),
              'components': {}, 'nonadditive_component_medians': True}
    for field in COMPONENTS:
        vals = [r[field] for r in rows if r[field] is not None]
        complete = len(vals) == len(rows) and bool(rows)
        result['components'][field] = {'status': 'available' if complete else 'unavailable',
            'available': len(vals), 'median': statistics.median(vals) if complete else None,
            'min': min(vals) if complete else None, 'max': max(vals) if complete else None,
            'range': max(vals) - min(vals) if complete else None,
            'missing_members': [r['case_id'] for r in rows if r[field] is None]}
    for method, field in [('native', COMPONENTS[0]), ('composite', COMPONENTS[-1])]:
        result[method + '_decision'] = decision(result['components'][field]['median'], frozen[method])
    return result


def counts(rows, method):
    calls = Counter(r[method + '_decision'] for r in rows)
    correct = sum(r[method + '_decision'] == r['expected_class'] + '-supported' for r in rows)
    wrong = sum(r[method + '_decision'] in ('Ca-supported', 'La-supported') and
                r[method + '_decision'] != r['expected_class'] + '-supported' for r in rows)
    return {'denominator': len(rows), 'correct': correct, 'wrong': wrong,
            'inconclusive': calls['inconclusive'], 'unavailable': calls['unavailable'], 'calls': dict(calls)}


def state_audit(prepared):
    if prepared is None or prepared['status'] != 'prepared':
        return {'status': 'unavailable', 'reason': None if prepared is None else prepared.get('reason', prepared['status'])}
    core = read_json(verify(prepared['core']['parent']))
    context = read_json(verify(prepared['representations']['context']['preparation']))
    eps = prepared['core']['endpoints']; ceps = prepared['representations']['context']['endpoints']
    elements = [r[0] for r in xyz(verify(eps['Ca']['xyz']))][1:]
    context_elements = [r[0] for r in xyz(verify(ceps['Ca']['xyz']))][1:]
    fragments = [{k: f.get(k) for k in ('kind', 'role', 'id', 'formal_charge', 'microstate_id', 'atom_count')}
                 for f in core['qm_fragments']]
    core_state = {'composition_without_metal': dict(Counter(elements)),
                  'ordered_elements_without_metal': elements,
                  'charges': {z: eps[z]['charge'] for z in ('Ca', 'La')},
                  'multiplicities': {z: eps[z]['multiplicity'] for z in ('Ca', 'La')},
                  'fragments': fragments, 'water_policy': core['fixed_core']['water_policy']}
    added = [{k: f.get(k) for k in ('kind', 'source', 'formal_charge')} for f in context['added_fragments']]
    return {'status': 'available', 'core_state_key': cache_key(core_state), 'core_state': core_state,
            'context_membership_key': cache_key(added), 'context_added_fragments': added,
            'context_added_formal_charge': context['added_formal_charge'],
            'context_atoms': len(context_elements) + 1,
            'context_composition_without_metal': dict(Counter(context_elements)),
            'context_charges': {z: ceps[z]['charge'] for z in ('Ca', 'La')},
            'core_source_coordinates_unchanged': context['core_source_coordinates_unchanged'],
            'water_inventory_unchanged': context['water_inventory_unchanged'],
            'core_preparation': prepared['core']['parent'],
            'context_preparation': prepared['representations']['context']['preparation']}


def summarize(sources, rows, frozen):
    groups = source_groups(sources); summaries = []
    for rep in REPS:
        index = {r['case_id']: r for r in rows if r['representation'] == rep}
        for root, members in groups.items():
            arms = {'La4': [m for m in members if m['source_conditioning_metal'] == 'La' and not m['canonical_coordinate_match']],
                    'Ca5': [m for m in members if m['source_conditioning_metal'] == 'Ca']}
            pools = {name: pool_summary(ms, index, frozen[rep]) for name, ms in arms.items()}
            balanced = {}; conditioning = {}
            for field in COMPONENTS:
                la = pools['La4']['components'][field]['median']; ca = pools['Ca5']['components'][field]['median']
                available = la is not None and ca is not None
                balanced[field] = (la + ca) / 2 if available else None
                conditioning[field] = ca - la if available else None
            item = {'root_case_id': root, 'representation': rep, 'expected_class': members[0]['expected_class'],
                    'biological_group': members[0]['biological_group'], 'pools': pools,
                    'equal_mean_of_arm_medians': balanced,
                    'Ca5_median_minus_La4_median': conditioning,
                    'canonical_replay_case_id': next(m['case_id'] for m in members if m['canonical_coordinate_match'])}
            for method, field in [('native', COMPONENTS[0]), ('composite', COMPONENTS[-1])]:
                item[method + '_decision'] = decision(balanced[field], frozen[rep][method])
            audited = [index[m['case_id']]['state_audit'] for m in members]
            available_states = [a for a in audited if a['status'] == 'available']
            item['state_accounting'] = {
                'prepared_sources': len(available_states), 'declared_sources': len(members),
                'core_state_invariant': len({a['core_state_key'] for a in available_states}) == 1 if len(available_states) == len(members) else None,
                'distinct_core_states': len({a['core_state_key'] for a in available_states}),
                'distinct_context_memberships': len({a['context_membership_key'] for a in available_states}),
                'context_added_charges': sorted({a['context_added_formal_charge'] for a in available_states}),
                'context_atom_counts': sorted({a['context_atoms'] for a in available_states}),
                'context_geometry_only_interpretation': False}
            summaries.append(item)
    totals = {}
    for rep in REPS:
        rr = [r for r in rows if r['representation'] == rep]
        subsets = {'all250': rr, 'canonical25_replay': [r for r in rr if r['canonical_coordinate_match']],
                   'primary225': [r for r in rr if r['primary_evaluation_pool']],
                   'La100_noncanonical': [r for r in rr if r['source_conditioning_metal'] == 'La' and not r['canonical_coordinate_match']],
                   'Ca125': [r for r in rr if r['source_conditioning_metal'] == 'Ca']}
        totals[rep] = {'single_sources': {name: {m: counts(rs, m) for m in ('native', 'composite')} for name, rs in subsets.items()},
                       'ensemble_descriptors': {}}
        parents = [s for s in summaries if s['representation'] == rep]
        for name in ('La4', 'Ca5', 'balanced'):
            selected = parents if name == 'balanced' else [dict(p['pools'][name], expected_class=p['expected_class']) for p in parents]
            totals[rep]['ensemble_descriptors'][name] = {m: counts(selected, m) for m in ('native', 'composite')}
    return summaries, totals


def compare(sources, preparation, inventory, collections, reference, output):
    source = read_json(sources); source_groups(source)
    prep = read_json(preparation); inv = read_json(inventory); ref = read_json(reference); frozen = bands(ref)
    if prep['source_manifest'] != record(sources) or inv['source_preparation'] != record(preparation):
        raise InvalidArtifact('source/preparation/inventory chain differs')
    inventory_pin = record(inventory); declared = {c['case_id']: c for c in source['cases']}
    reference_rows = {(r['case_id'], r['representation']): r for r in ref['rows']}
    available = {c['case_id']: c for c in inv['cases']}; prepared = {c['case_id']: c for c in prep['cases']}
    if not set(available) <= set(declared) or not set(prepared) <= set(declared):
        raise InvalidArtifact('undeclared source in prepared inventory')
    for c in source['cases']:
        original = reference_rows[(c['root_case_id'], 'core')]
        if any(c[k] != original[k] for k in ('expected_class', 'biological_group')):
            raise InvalidArtifact('source label/group differs from the frozen parent reference')
        if c['case_id'] in available and any(available[c['case_id']][k] != c[k] for k in
                ('root_case_id', 'source_conditioning_metal', 'canonical_coordinate_match', 'expected_class')):
            raise InvalidArtifact('native inventory source identity differs')
    lookup = {}; failures = []; collection_pins = []; reuse_audits = []
    for path in collections:
        col = read_json(path); collection_pins.append(record(path))
        if col['inventory'] != inventory_pin or col['protocol_id'] != COMPOSITE_PROTOCOL:
            raise InvalidArtifact('collection inventory/method differs')
        for row in col['rows']:
            if row['solver'] != 'native': continue
            key = (row['case_id'], row['representation'], row['metal'])
            if key[0] not in available or key[1] not in REPS or key[2] not in ('Ca', 'La'):
                raise InvalidArtifact('undeclared transfer endpoint')
            endpoint = available[key[0]]['representations'][key[1]]['endpoints'][key[2]]
            reused = checked_scoped_transfer(row, endpoint, inventory_pin, col['manifest'])
            if reused:
                reuse_audits.append({'case_id': key[0], 'representation': key[1], 'metal': key[2], 'endpoints': reused})
            if row['status'] != 'complete':
                failures.append({'collection': record(path), 'row': row}); lookup.setdefault(key, row); continue
            if key in lookup and lookup[key]['status'] == 'complete' and row != lookup[key]:
                raise InvalidArtifact('conflicting duplicate transfer result')
            lookup[key] = row
    rows = []
    for case in source['cases']:
        cid = case['case_id']; old = available.get(cid); audit = state_audit(prepared.get(cid))
        for rep in REPS:
            row = {k: case[k] for k in ('case_id', 'root_case_id', 'source_conditioning_metal', 'canonical_coordinate_match',
                    'primary_evaluation_pool', 'biological_group', 'expected_class', 'evidence_stratum')}
            row.update(representation=rep, state_audit=audit, source_structure=case['source_structure'],
                       preparation_status=prepared.get(cid, {}).get('status', 'missing'),
                       preparation_reason=prepared.get(cid, {}).get('reason'), native_status='unavailable', composite_status='unavailable',
                       low_level_endpoints={}, **{field: None for field in COMPONENTS})
            if old is not None and rep in old['representations']:
                eps = old['representations'][rep]['endpoints']; row['native_MACE_endpoints'] = eps
                current = prepared[cid]['core']['endpoints'] if rep == 'core' else prepared[cid]['representations'][rep]['endpoints']
                row['scored_vs_current_preparation_max_displacement_A'] = {}
                for z, ep in eps.items():
                    actual = native_endpoint(read_json(verify(ep['native_MACE_receipt'])), inv['model'])
                    if actual != ep: raise InvalidArtifact('native endpoint differs from original receipt/state')
                    aa = xyz(verify(ep['xyz'])); bb = xyz(verify(current[z]['xyz']))
                    if (len(aa) != len(bb) or [a[0] for a in aa] != [b[0] for b in bb]
                            or ep['charge'] != current[z]['charge'] or ep['multiplicity'] != current[z]['multiplicity']):
                        raise InvalidArtifact('native endpoint differs from the prepared source physical state')
                    displacement = max(abs(a[i] - b[i]) for a, b in zip(aa, bb) for i in (1, 2, 3))
                    tolerance = prep.get('numeric_reconciliation', {}).get('tolerance_A', 0) if case['canonical_coordinate_match'] else 0
                    if displacement > tolerance:
                        raise InvalidArtifact('native endpoint geometry differs from the prepared source')
                    row['scored_vs_current_preparation_max_displacement_A'][z] = displacement
                row['native_R_model_kcal_mol'] = old['representations'][rep]['native_R_model_kcal_mol']
                row['native_status'] = 'available'
                ll = {z: lookup.get((cid, rep, z)) for z in ('Ca', 'La')}; row['low_level_endpoints'] = ll
                for medium, field in [('vacuum', 'GFN2_vacuum_R_kcal_mol'), ('alpb', 'GFN2_ALPB_R_kcal_mol')]:
                    if all(v is not None and v['endpoints'][medium]['status'] == 'complete' for v in ll.values()):
                        row[field] = (ll['Ca'][medium + '_hartree'] - ll['La'][medium + '_hartree']) * HA_TO_KCAL
                row['endpoint_solvation_kcal_mol'] = {
                    z: ll[z]['delta_solv_hartree'] * HA_TO_KCAL if ll[z] is not None and ll[z]['status'] == 'complete' else None
                    for z in ('Ca', 'La')}
                if all(v is not None and v['status'] == 'complete' for v in ll.values()):
                    result = mix_pair({z: eps[z]['native_MACE_energy_eV'] for z in ('Ca', 'La')},
                                      {z: ll[z]['vacuum_hartree'] for z in ('Ca', 'La')},
                                      {z: ll[z]['alpb_hartree'] for z in ('Ca', 'La')})
                    if result['native_R_model_kcal_mol'] != row['native_R_model_kcal_mol']:
                        raise InvalidArtifact('native contrast algebra differs')
                    row.update(result, composite_status='available')
            for method, field in [('native', COMPONENTS[0]), ('composite', COMPONENTS[-1])]:
                row[method + '_decision'] = decision(row[field], frozen[rep][method])
            if case['canonical_coordinate_match']:
                original = reference_rows[(case['root_case_id'], rep)]
                row['canonical_replay_comparison'] = {
                    'original_case_id': case['root_case_id'],
                    'differences_model_kcal_mol': {field: row[field] - original[field] if row[field] is not None else None
                                                  for field in COMPONENTS},
                    'original_composite_decision': original['composite_decision'],
                    'current_composite_decision': row['composite_decision']}
            rows.append(row)
    summaries, totals = summarize(source, rows, frozen)
    result = {'schema_version': SCHEMA, 'sources': record(sources), 'preparation': record(preparation),
              'inventory': inventory_pin, 'collections': collection_pins, 'frozen_reference': record(reference),
              'frozen_bands': frozen, 'rows': rows, 'protein_summaries': summaries, 'counts': totals,
              'source_denominator': 250, 'protein_denominator': 25, 'case_representation_denominator': 500,
              'failed_transfer_history': failures, 'native_inventory_failures': inv.get('failures', []),
              'explicit_original_receipt_reuse': reuse_audits,
              'threshold_refitted': False, 'old_bands_used_as_transfer_test_only': True,
              'ensemble_interpretation': 'declared structural descriptor; not a thermal population or free energy',
              'new_scientific_calls': 0, 'production_changed': False, 'all_evidence_consumed': True,
              'implementation': record(__file__)}
    write_new(output, result)
    return result


def report(result, output):
    r = read_json(result)
    lines = ['# Frozen fold-conditioning robustness comparison', '',
             'Existing 25 PQQ proteins, 250 saved source geometries. All evidence was consumed; structural samples are not independent biological labels.', '',
             '| Representation | Method | Descriptor | Correct | Wrong | Inconclusive | Unavailable | Denominator |',
             '|---|---|---|---:|---:|---:|---:|---:|']
    for rep, tt in r['counts'].items():
        for kind, groups in tt.items():
            for name, methods in groups.items():
                for method, c in methods.items():
                    lines.append(f"| {rep} | {method} | {name} | {c['correct']} | {c['wrong']} | {c['inconclusive']} | {c['unavailable']} | {c['denominator']} |")
    primary_rows = [x for x in r['rows'] if x['representation'] == 'core']
    failures = Counter(x['preparation_reason'] or x['preparation_status'] for x in primary_rows
                       if x['preparation_status'] != 'prepared')
    states = [s['state_accounting'] for s in r['protein_summaries'] if s['representation'] == 'context']
    lines += ['', '## State accounting', '',
              f"Prepared sources: {sum(x['preparation_status'] == 'prepared' for x in primary_rows)}/250.",
              f"Preparation failures: {dict(failures)}.",
              f"Complete protein pools with invariant core state: {sum(s['core_state_invariant'] is True for s in states)}; "
              f"incomplete pools: {sum(s['core_state_invariant'] is None for s in states)}; "
              f"complete pools with differing core states: {sum(s['core_state_invariant'] is False for s in states)}.",
              f"Proteins with multiple observed context memberships: {sum(s['distinct_context_memberships'] > 1 for s in states)}/25."]
    lines += ['', 'Bands are frozen from the existing representation-specific calibration. No thresholds were fit. The original 25 canonical source replays are separate from the primary 225 source transfer geometries.', '',
              'La4/Ca5 are medians of all four noncanonical La-conditioned/all five Ca-conditioned sources per protein. Balanced is the equal mean of those two medians. Missing members make the corresponding descriptor unavailable; no successful-member selection occurs.', '',
              'Native OMOL, GFN2 vacuum/ALPB contrasts and transfer corrections remain separate in the JSON. Component medians do not add to the median composite in general. Context membership can change with geometry; fixed-core state invariants and context composition changes are explicitly recorded.', '',
              'This report is an existing-band transfer test, not a newly calibrated ensemble classifier, broad biological validation, thermal population or binding free energy.', '',
              f"Actual collection: `{Path(result).resolve()}`", '']
    path = Path(output)
    if path.exists(): raise InvalidArtifact('immutable report already exists')
    path.write_text('\n'.join(lines))


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='op', required=True)
    q = sub.add_parser('compare')
    for name in ('sources', 'preparation', 'inventory', 'reference', 'output'): q.add_argument('--' + name, required=True)
    q.add_argument('--collection', action='append', default=[], dest='collections')
    q = sub.add_parser('report'); q.add_argument('--result', required=True); q.add_argument('--output', required=True)
    a = vars(p.parse_args()); op = a.pop('op'); result = globals()[op](**a)
    if op == 'compare': print(json.dumps({'output': a['output'], 'counts': result['counts']}))


if __name__ == '__main__': main()
