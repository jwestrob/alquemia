"""Overlay two actually completed origin-based searches on the full225 ledger."""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path

from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from adaptive_minimal_pool import PROTOCOL as BASE_PROTOCOL, checked_reference
from adaptive_origin_recovery import CASES, CANDIDATES, PROTOCOL, validate_pool
from nikasha_pool import choose_rows
from nikasha_pool_compare import aggregate, method_record
from accommodation_fold_proposals import summarize_rows


def compare(base, recovery, output):
    original = read_json(base); fresh = read_json(recovery)
    m = read_json(verify(fresh['manifest'])); validate_pool(verify(fresh['manifest']))
    ref = read_json(verify(original['reference']))
    checked_reference(verify(original['reference']), {'signature': ref['signature']})
    if (original['protocol_id'] != BASE_PROTOCOL or fresh['protocol_id'] != PROTOCOL
            or m['reference'] != original['reference'] or len(original['rows']) != 225
            or len(fresh['cases']) != 2 or {r['case_id'] for r in fresh['cases']} != set(CASES)):
        raise InvalidArtifact('recovery/source/reference scope mismatch')
    updates = {r['case_id']: r for r in fresh['cases']}
    bands = {**original['bands'], 'minimal_recovered': ref['variants']['operational']['bands']}
    rows = []
    for old in original['rows']:
        row = {**old, 'methods': dict(old['methods']), 'recovery_result': None}
        value = old['methods']['minimal_operational']['R']
        if old['case_id'] in updates:
            new = updates[old['case_id']]
            if value is not None:
                raise InvalidArtifact('recovery would replace a previously available score')
            if new['status'] == 'prepared':
                if [r['id'] for r in new['candidates']] != list(CANDIDATES):
                    raise InvalidArtifact('not the declared three-candidate pool')
                if choose_rows(new['matrix'], CANDIDATES) != new['pool']:
                    raise InvalidArtifact('saved new energy/selection algebra differs')
            elif new['pool']['status'] != 'unavailable':
                raise InvalidArtifact('failed proposal was silently scored')
            value = new['pool']['operational']['composite_R_model_kcal_mol'] if new['pool']['operational'] else None
            row['recovery_result'] = {'collection': record(recovery), 'case': new}
        row['methods']['minimal_recovered'] = method_record(value, bands['minimal_recovered'], old['expected_class'])
        rows.append(row)
    index = {r['case_id']: r for r in rows}; pools = []; triples = []
    for group in original['pools']:
        if group['pool'] == 'balanced':
            values = [next(g for g in pools if g['root_case_id'] == group['root_case_id'] and g['pool'] == name)
                      ['methods']['minimal_recovered']['R'] for name in ('La4', 'Ca5')]
            value = sum(values)/2 if all(v is not None for v in values) else None
            new = {**method_record(value, bands['minimal_recovered'], group['expected_class']),
                   'Ca5_minus_La4': values[1]-values[0] if value is not None else None}
        else:
            new = aggregate(group['members'], index, 'minimal_recovered', bands['minimal_recovered'], group['expected_class'])
        pools.append({**group, 'methods': {**group['methods'], 'minimal_recovered': new}})
    for group in original['triples']:
        new = aggregate(group['members'], index, 'minimal_recovered', bands['minimal_recovered'], group['expected_class'])
        triples.append({**group, 'methods': {**group['methods'], 'minimal_recovered': new}})
    subsets = {'all225': rows, 'La100': [r for r in rows if r['source_conditioning_metal'] == 'La'],
               'Ca125': [r for r in rows if r['source_conditioning_metal'] == 'Ca'], 'all100_triples': triples,
               **{p: [r for r in pools if r['pool'] == p] for p in ('La4', 'Ca5', 'balanced')}}
    counts = {name: summarize_rows(values, bands) for name, values in subsets.items()}
    matched = {}
    for name, values in subsets.items():
        matched[name] = {}
        for other in ('released', 'full_adaptive', 'minimal_operational'):
            common = [r for r in values if all(r['methods'][k]['decision'] != 'unavailable'
                      for k in (other, 'minimal_recovered'))]
            matched[name][other] = {'common': len(common),
                'counts': summarize_rows(common, (other, 'minimal_recovered')),
                'transitions': dict(Counter(r['methods'][other]['outcome'] + '->' +
                    r['methods']['minimal_recovered']['outcome'] for r in common))}
    result = {'protocol_id': PROTOCOL, 'base_comparison': record(base), 'recovery': record(recovery),
        'reference': original['reference'], 'implementation': record(__file__),
        'bands': bands, 'rows': rows, 'pools': pools, 'triples': triples, 'counts': counts,
        'matched': matched, 'new_molecular_calls_in_comparison': 0,
        'changed_source_ids': list(CASES), 'thresholds_refitted': False,
        'historical_unavailable_rows_preserved': True, 'production_changed': False}
    write_new(output, result)
    return {'all225': counts['all225']['minimal_recovered'], 'matched_released': matched['all225']['released']}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('base', 'recovery', 'output'):
        p.add_argument('--' + name, required=True, type=Path)
    print(json.dumps(compare(**vars(p.parse_args())), indent=2))


if __name__ == '__main__':
    main()
