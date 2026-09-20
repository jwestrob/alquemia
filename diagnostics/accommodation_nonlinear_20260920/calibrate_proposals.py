"""Canonical-only developmental reference for the declared proposal descriptor."""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from accommodation_proposals import PROTOCOL, SETTINGS, score, select_endpoint
from compact_solvation_compare import MINIMUM_CALIBRATION_GAP
from compact_solvation_scanner import decision

REFERENCE_ID = 'PQQ_MACE_donor_proposal_composite_canonical_reference_v1'


def build(result, manifest, original):
    if result['protocol_id'] != PROTOCOL or manifest['protocol_id'] != PROTOCOL:
        raise InvalidArtifact('proposal method differs')
    if manifest['settings'] != SETTINGS:
        raise InvalidArtifact('geometry-selection policy differs')
    # Membership and labels come from the original designated calibration record,
    # never from an advantageous score or a prediction on an unknown PLM target.
    canonical = [r for r in original['rows']
                 if r['representation'] == 'context' and r['role'] == 'calibration']
    ids = {r['case_id'] for r in canonical}
    if len(canonical) != 25 or len(ids) != 25:
        raise InvalidArtifact('canonical calibration membership differs')
    cases = {c['case_id']: c for c in result['cases']}
    if len(cases) != len(result['cases']) or not ids <= set(cases):
        raise InvalidArtifact('duplicate or missing calibration identity')
    endpoints = {(r['case_id'], r['metal']): r for r in result['endpoints']}
    if len(endpoints) != len(result['endpoints']):
        raise InvalidArtifact('duplicate endpoint identity')
    tasks = {(r['case_id'], r['metal']): r for r in manifest['tasks']}
    rows = []
    for old in canonical:
        c = cases[old['case_id']]
        if c['expected_class'] != old['expected_class'] or c['label_scope'] != old['label_scope']:
            raise InvalidArtifact('canonical biological class or evidence stratum differs')
        row = {'case_id': c['case_id'], 'expected_class': old['expected_class'],
               'biological_group': old['biological_group'], 'status': 'unavailable',
               'R_model_kcal_mol': None, 'source_endpoints': {}, 'selection': {}}
        if c['status'] == 'available':
            pair = {}
            for metal in ('Ca', 'La'):
                e = endpoints[c['case_id'], metal]
                if e['status'] != 'available':
                    raise InvalidArtifact('case marked available with missing endpoint')
                sel = select_endpoint(e['origin_components'], e['proposal_components'])
                if sel != e['selection']:
                    raise InvalidArtifact('actual endpoint selection differs')
                pair[metal] = sel['selected_components']
                task = tasks[c['case_id'], metal]
                row['source_endpoints'][metal] = {k: task[k] for k in
                    ('xyz', 'mapping', 'source_preparation', 'charge', 'multiplicity')}
                row['selection'][metal] = {'selected': sel['selected'],
                    'xyz': e['proposal_xyz'] if sel['selected'] == 'proposal' else e['origin_xyz'],
                    'proposal_receipt': e['proposal_receipt']}
            computed = score(pair['Ca'], pair['La'])
            if computed != c['R_selected']:
                raise InvalidArtifact('collected score differs from actual component algebra')
            row.update(status='available', R_model_kcal_mol=computed['composite_R_model_kcal_mol'])
        rows.append(row)
    available = [r for r in rows if r['status'] == 'available']
    bands = None; gap = None; extrema = None
    status = 'unavailable_calibration_member'
    if len(available) == 25:
        groups = {z: [r['R_model_kcal_mol'] for r in rows if r['expected_class'] == z]
                  for z in ('Ca', 'La')}
        if not all(groups.values()):
            raise InvalidArtifact('missing calibration class')
        extrema = {'Ca_max': max(groups['Ca']), 'La_min': min(groups['La'])}
        gap = extrema['La_min'] - extrema['Ca_max']
        status = 'available' if gap > MINIMUM_CALIBRATION_GAP else 'unsupported_class_separation'
        if status == 'available':
            bands = {'Ca_supported_max_R_model_kcal_mol': extrema['Ca_max'],
                     'La_supported_min_R_model_kcal_mol': extrema['La_min'],
                     'protocol_id': PROTOCOL, 'reference_id': REFERENCE_ID,
                     'representation': 'context', 'solver': 'native'}
    transferred = []
    for c in result['cases']:
        value = c['R_selected']['composite_R_model_kcal_mol'] if c['R_selected'] else None
        transferred.append({'case_id': c['case_id'],
            'role': 'calibration' if c['case_id'] in ids else
                    'unlabeled_prediction' if c['expected_class'] is None else 'consumed_crystal_transfer',
            'expected_class': c['expected_class'], 'R_model_kcal_mol': value,
            'old_band_decision': c['old_band_selected_decision'],
            'new_developmental_reference_decision': decision(value, bands) if bands else None})
    return {'reference_id': REFERENCE_ID, 'protocol_id': PROTOCOL, 'status': status,
            'calibration_denominator': 25, 'available_calibration': len(available),
            'calibration_class_counts': dict(Counter(r['expected_class'] for r in rows)),
            'rows': rows, 'bands': bands, 'class_extrema': extrema,
            'gap_model_kcal_mol': gap, 'minimum_gap_model_kcal_mol': MINIMUM_CALIBRATION_GAP,
            'rule': 'identical class-extrema rule and minimum gap as prior composite calibration',
            'settings': manifest['settings'], 'model': manifest['model'],
            'orca': manifest['orca'], 'source_implementation': manifest['implementation'],
            'old_frozen_bands': manifest['frozen_context_bands'],
            'initial_source_decisions': transferred, 'production_changed': False,
            'noncanonical_folds_used_for_calibration': False,
            'crystals_or_PLM_used_for_calibration': False,
            'new_molecular_calls': 0,
            'interpretation': 'new developmental calibration; calibration success is not validation'}


def prepare(result, agreement, output):
    d = read_json(result); mp = verify(d['manifest']); m = read_json(mp)
    oldp = verify(m['frozen_comparison']); old = read_json(oldp)
    for c in m['cases']:
        for pin in c['maps'].values(): verify(pin)
        for e in c['origins'].values(): verify(e['xyz'])
    r = build(d, m, old)
    r.update(result=record(result), manifest=record(mp), original_reference=record(oldp),
             agreement=record(agreement), implementation=record(__file__),
             original_rule_implementation=record(ROOT / 'scripts/compact_solvation_compare.py'))
    write_new(output, r)
    return {k: r[k] for k in ('reference_id', 'status', 'bands', 'gap_model_kcal_mol',
                              'calibration_denominator', 'available_calibration')}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for k in ('result', 'agreement', 'output'): p.add_argument('--' + k, required=True)
    print(json.dumps(prepare(**vars(p.parse_args())), indent=2))
