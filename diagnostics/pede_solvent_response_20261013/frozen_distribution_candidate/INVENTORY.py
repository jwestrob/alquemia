#!/usr/bin/env python3
"""Metadata inventory only: no coordinates, populations, energies or scores evaluated."""
import hashlib
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_inventory_v1'
ORDER = ['origin', 'adaptive_Ca', 'adaptive_La']
SOURCES = {
    'expanded135': 'workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3',
    'primary_ca125': 'workspaces/orca_parameter_audit_20261012/ca_conditioned_gfn2_slurm_v1',
    'mxaf_recovery2': 'workspaces/ca_coverage_recovery_20261012/mxaf_solvent_v1',
    'disulfide_recovery1': 'workspaces/ca_coverage_recovery_20261012/disulfide/solvent_v1',
    'two_ca_q9_5': 'workspaces/ca_coverage_recovery_20261012/two_ca_solvent_q9_v1',
    'two_ca_c5q88_10': 'workspaces/ca_coverage_recovery_20261012/two_ca_solvent_c5q88_v1',
}


@lru_cache(None)
def read(path):
    return json.loads(Path(path).read_text())


@lru_cache(None)
def pin(path):
    path = Path(path).resolve()
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return {'path': str(path), 'sha256': h.hexdigest()}


def checked(p):
    assert pin(p['path']) == p, ('changed artifact', p['path'])
    return read(p['path'])


def ident(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def native_endpoint(component, medium):
    pp = component['inputs'][medium]
    point = checked(pp)
    backend = point['backends']['xtb']
    evidence = component['verified_output_evidence'][medium]['xtb']
    assert evidence['execution'] == backend['execution']
    receipt = checked(evidence['execution'])
    assert backend['status'] == 'available' and receipt['returncode'] == 0
    outputs = {Path(p['path']).name: p for p in evidence['outputs']}
    names = {'charges': 'charges', 'stdout': 'stdout.txt', 'executed_xyz': 'core.xyz'}
    result = {'point': pp, 'execution': evidence['execution'], 'task': point['task'],
              'point_manifest': point['manifest'], 'status': backend['status']}
    for field, name in names.items():
        p = outputs[name]
        assert pin(p['path']) == p
        result[field] = p
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    assert not (OUT / 'INVENTORY.json').exists(), 'frozen inventory already exists'
    records, tasks, native_keys, source_pins, case_lookup = [], [], {}, {}, {}
    triples = []
    for family, directory in SOURCES.items():
        base = ROOT / directory
        mp, cp = pin(str(base / 'manifest.json')), pin(str(base / 'FINAL_COLLECTION.json'))
        m, collection = checked(mp), checked(cp)
        assert collection['manifest'] == mp
        source_pins[family] = {'manifest': mp, 'collection': cp, 'runtime': m['runtime'],
                               'model': m['model'], 'settings': m['settings']}
        outcomes = {(r['source_case_id'], r['case_id']): (i, r) for i, r in enumerate(collection['rows'])}
        for index, c in enumerate(m['cases']):
            c = {'case_id': None} | c
            sid = c['source_case_id']
            ri, old = outcomes[sid, c['case_id']]
            assert old['case_id'] == c['case_id']
            if c.get('execution_status') == 'reuse_actual_corrected_pool':
                prior = checked(c['prior_collection'])
                prior_manifest = checked(prior['manifest'])
                same = [a for a in prior_manifest['cases'] if a['case_id'] == c['case_id']]
                assert len(same) == 1
                assert all(c[k] == same[0][k] for k in ('state_signature', 'maps', 'source_preparation'))
                c = c | {'matrix': same[0]['matrix']}
            rid = family + '::' + sid + '::' + str(c['selection_id'])
            r = {k: c[k] for k in ('case_id', 'source_case_id', 'selection_id', 'protein',
                  'biological_group', 'expected_class', 'source_conditioning_metal')}
            r.update(record_id=rid, family=family,
                     role=c.get('role', family), original_source_id=c.get('original_source_id', sid),
                     baseline_manifest=mp, baseline_manifest_case_index=index,
                     baseline_collection=cp, baseline_collection_row_index=ri,
                     baseline_pool_status=old['pool']['status'],
                     preparation_status=c.get('preparation_status', 'archived_admitted'),
                     preparation_reason=c.get('preparation_reason'),
                     matrix={'Ca': {}, 'La': {}}, missing_cells=[])
            for key in ('state_signature', 'maps', 'source_preparation', 'prepared_admission',
                        'native_source_record', 'new_source_admission', 'original_exclusion'):
                if key in c:
                    r[key] = c[key]
            for metal in ('Ca', 'La'):
                for candidate in ORDER:
                    cell = c.get('matrix', {}).get(metal, {}).get(candidate)
                    component = old.get('components', {}).get(metal + '__' + candidate)
                    reason = None
                    if not cell:
                        reason = 'missing_native_cell_or_preparation'
                    elif not component or not component.get('accepted'):
                        reason = 'unavailable_or_unqualified_archived_solvent_pair'
                    if reason:
                        missing = {'metal': metal, 'candidate': candidate, 'reason': reason}
                        r['missing_cells'].append(missing)
                        r['matrix'][metal][candidate] = {'status': 'unavailable', 'reason': reason}
                        continue
                    endpoints = {medium: native_endpoint(component, medium) for medium in ('vacuum', 'water')}
                    a, b = endpoints['vacuum']['task'], endpoints['water']['task']
                    fields = ('metal', 'charge', 'multiplicity', 'valence_electrons', 'atom_count', 'xyz')
                    assert all(a[k] == b[k] for k in fields)
                    assert a['metal'] == metal and a['multiplicity'] == 1
                    assert a['xyz'] == cell['xyz'] and a['charge'] == cell['charge']
                    mace = checked(cell['MACE'])
                    assert mace['energy_eV'] == cell['MACE_eV']
                    assert pin(cell['xyz']['path']) == cell['xyz']
                    # Exact-byte dedup only; no coordinates, populations or energies recomputed.
                    scientific_key = {k: a[k] for k in fields if k != 'xyz'} | {
                        'xyz_sha256': a['xyz']['sha256'],
                        'state_signature_sha256': c['state_signature']['sha256'],
                        'vacuum_charges_sha256': endpoints['vacuum']['charges']['sha256'],
                        'water_charges_sha256': endpoints['water']['charges']['sha256'],
                        'runtime': m['runtime'], 'settings': m['settings']}
                    key = ident(scientific_key)
                    tid = 'native_' + key
                    if key not in native_keys:
                        t = {'task_id': tid, 'scientific_input_key': scientific_key,
                             **{k: a[k] for k in fields},
                             'state_signature': c['state_signature'], 'maps': c['maps'],
                             'source_preparation': c['source_preparation'], 'endpoints': endpoints,
                             'archived_water_components': {
                                 'status': 'read_from_pinned_stdout_during_Slurm_preflight',
                                 'stdout': endpoints['water']['stdout'],
                                 'required_terms': ['Gelec', 'Ghb', 'Gsasa', 'Gshift']},
                             'uses': []}
                        native_keys[key] = t
                        tasks.append(t)
                    use = {'record_id': rid, 'source_case_id': sid, 'case_id': c['case_id'],
                           'candidate': candidate, 'metal': metal, 'MACE': cell['MACE'],
                           'MACE_eV': cell['MACE_eV']}
                    native_keys[key]['uses'].append(use)
                    r['matrix'][metal][candidate] = {'status': 'available_archived_inputs',
                        'native_task_id': tid, 'xyz': cell['xyz'], 'MACE': cell['MACE'],
                        'MACE_eV': cell['MACE_eV']}
            r['input_status'] = 'complete_pool' if not r['missing_cells'] else 'incomplete_pool'
            r['origin_input_status'] = ('available' if all(r['matrix'][z]['origin']['status'] ==
                'available_archived_inputs' for z in ('Ca', 'La')) else 'unavailable')
            records.append(r)
            case_lookup[family, sid] = r
        if family == 'expanded135':
            triples = m['triples']

    merged_pin = pin(str(ROOT / 'workspaces/ca_coverage_recovery_20261012/two_ca_scoring_final_v1/RESULT.json'))
    merged = checked(merged_pin)
    current = []
    recovery_family = {'technical_recovery': 'mxaf_recovery2',
                       'preparation_recovery': 'disulfide_recovery1'}
    for index, row in enumerate(merged['rows']):
        primary = row['primary']
        active = [(k, v) for k, v in row.items() if k != 'primary' and v is not None]
        assert len(active) <= 1
        family, selected = 'primary_ca125', primary
        kind = 'unchanged_primary'
        if active:
            kind, selected = active[0]
            family = ('two_ca_q9_5' if selected['biological_group'] == 'q9z4j7' else
                      'two_ca_c5q88_10') if kind == 'two_Ca_source_recovery' else recovery_family[kind]
        record = case_lookup[family, selected['source_case_id']]
        assert record['case_id'] == selected['case_id']
        current.append({'original_source_id': primary['source_case_id'], 'record_id': record['record_id'],
                        'selection_kind': kind, 'actual_source_case_id': record['source_case_id'],
                        'baseline_overlay': merged_pin, 'baseline_overlay_row_index': index})
    assert len(current) == 125
    canonical = [r['record_id'] for r in records if r['family'] == 'expanded135' and r['role'] == 'canonical_calibration']
    assert len(canonical) == 25
    canonical_members = set(canonical)
    for t in tasks:
        t['stage'] = 'canonical' if any(u['record_id'] in canonical_members for u in t['uses']) else 'transfer'
    pools = {r['record_id']: r for r in records}
    summary = {'records_by_family': dict(Counter(r['family'] for r in records)),
               'complete_pools_by_family': dict(Counter(r['family'] for r in records if r['input_status'] == 'complete_pool')),
               'unique_native_inputs_by_stage': dict(Counter(t['stage'] for t in tasks)),
               'native_inputs': len(tasks), 'declared_component_replays_per_native_input': 2,
               'state_pairs': sorted({(t['metal'], t['charge'], t['multiplicity']) for t in tasks}),
               'atom_counts': sorted({t['atom_count'] for t in tasks}),
               'current_overlay_complete': all(pools[r['record_id']]['input_status'] == 'complete_pool' for r in current),
               'scientific_evaluations': 0, 'candidate_scores_calculated': 0}
    m = {'protocol': 'frozen_distribution_candidate_inventory_v1',
         'status': 'metadata_ready_not_a_scientific_execution_manifest',
         'candidate_energy': 'E_OMOL,vac + G_nativeALPB(q_GFN2,vac)',
         'candidate_components': ['Gelec', 'Ghb', 'Gsasa', 'Gshift'],
         'inventory_implementation': pin(str(Path(__file__))), 'sources': source_pins,
         'baseline_canonical': pin(str(ROOT / SOURCES['expanded135'] / 'CANONICAL_COLLECTION.json')),
         'calibration_rule_source': pin(str(ROOT / 'diagnostics/orca_parameter_audit_20261012/expanded_transfer_inventory/REEVALUATE.py')),
         'calibration_rule': {'designated_record_ids': canonical, 'requires_all': 25,
             'minimum_gap_strictly_greater_kcal_mol': 0.02, 'Ca_max': 'max designated Ca scores',
             'La_min': 'min designated La scores', 'fit_roles': ['canonical_calibration'],
             'gate': 'operational canonical calibration available before dependent transfer'},
         'selection_rule': {'candidate_order': ORDER, 'origin_retained_unless_minimum_work_below_kcal_mol': -0.1,
             'mathematical_minimum_separate': True, 'full_pool_required': True,
             'origin_diagnostic_independently_available': True},
         'cohorts': {'expanded135': [r['record_id'] for r in records if r['family'] == 'expanded135'],
             'primary_ca125': [r['record_id'] for r in records if r['family'] == 'primary_ca125'],
             'current_ca125_annotated_overlay': current,
             'two_ca15': [r['record_id'] for r in records if r['family'].startswith('two_ca_')],
             'existing_100_triples': triples},
         'records': records, 'native_inputs': tasks, 'summary': summary,
         'limits': ['consumed development references; no fresh blind claim',
                   'retain primary failures and all six historical triple exclusions',
                   'no new Ca-conditioned triple/vote policy',
                   'no old-band inheritance or fitted response factor',
                   'no strict208/tenfold protocol mixing',
                   'no omitted unknown PLM labels inferred'],
         'production_changed': False}
    (OUT / 'INVENTORY.json').write_text(json.dumps(m, indent=2, allow_nan=False) + '\n')
    (OUT / 'SUMMARY.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({'inventory': pin(str(OUT / 'INVENTORY.json')), 'summary': summary}, indent=2))


if __name__ == '__main__':
    main()
