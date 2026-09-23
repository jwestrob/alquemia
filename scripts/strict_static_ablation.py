"""Separate static-pocket and accommodation utility using completed strict energies."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import json

from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from nikasha_pool import score
from nikasha_pool_compare import extrema_reference
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome, strict_summary, summarize_rows

PROTOCOL = 'Nikasha_tenfold_strict_static_ablation_v1'
METHOD = 'union_strict_static'
ADAPTIVE = 'union_precision_strict'


def origin_score(matrix):
    pair = [matrix[z]['origin'] for z in ('Ca', 'La')]
    if any(c['status'] != 'complete' or c['components'] is None for c in pair):
        return None
    return score(*(c['components'] for c in pair))['composite_R_model_kcal_mol']


def calibrate(qualification, output):
    q = read_json(qualification)
    if q['agreeing_cells'] != 384 or q['counts']['all32']['qualified'] != 32:
        raise InvalidArtifact('Completed uniform strict qualification required')
    rows = [{k: r[k] for k in ('case_id', 'role', 'expected_class')} |
            {'R_model_kcal_mol': origin_score(r['branches']['fresh']['matrix'])}
            for r in q['rows'] if r['role'] == 'calibration']
    if len({r['case_id'] for r in rows}) != 25:
        raise InvalidArtifact('Exactly25 distinct designated canonical inputs required')
    ref = extrema_reference(rows, 'origin', PROTOCOL)
    result = {'protocol_id': PROTOCOL, 'qualification': record(qualification),
              'adaptive_reference': q['reference'], 'static_reference': ref,
              'frozen_at_UTC': datetime.now(timezone.utc).isoformat(),
              'rule': 'Unchanged canonical25 extrema and minimum-gap rule, origin only.',
              'transfer_rows_used_for_fit': False, 'new_molecular_calls': 0,
              'implementation': record(__file__)}
    write_new(output, result)
    return {k: v for k, v in ref.items() if k != 'rows'}


def classify(value, bands, expected):
    call = decision(value, bands) if bands else 'unavailable'
    return {'R': value, 'decision': call, 'outcome': outcome(call, expected)}


def compare(transfer, reference, output):
    src, ref = read_json(transfer), read_json(reference)
    verify(ref['qualification'])
    if (ref['protocol_id'] != PROTOCOL or src['reference'] != ref['adaptive_reference']
            or src['qualification'] != ref['qualification']):
        raise InvalidArtifact('Strict fresh qualification/reference mismatch')
    band = ref['static_reference']['bands']
    if band is None:
        raise InvalidArtifact('Static canonical reference unavailable; no substituted bands')
    bands = {**src['bands'], METHOD: band}
    matrices = {c['case_id']: c['matrix'] for c in src['actual_pools']}
    if len(src['rows']) != 225 or len(matrices) != 208:
        raise InvalidArtifact('Original225/208 source denominator required')
    rows = []
    for old in src['rows']:
        value = origin_score(matrices[old['case_id']]) if old['case_id'] in matrices else None
        adapted = old['methods'][ADAPTIVE]['R']
        rows.append({**old, 'methods': {**old['methods'], METHOD:
                    classify(value, band, old['expected_class'])},
                     'strict_accommodation_delta_R':
                     adapted-value if adapted is not None and value is not None else None,
                     'strict_origin_under_adaptive_bands':
                     classify(value, bands[ADAPTIVE], old['expected_class'])})
    index = {r['case_id']: r for r in rows}
    groups = []
    for old in src['pools']:
        if old['pool'] in ('La4', 'Ca5'):
            value = strict_summary(old['members'], index, METHOD, band, old['expected_class'])
        else:
            arms = [next(p for p in groups if p['root_case_id'] == old['root_case_id']
                         and p['pool'] == name)['methods'][METHOD] for name in ('La4', 'Ca5')]
            val = None if any(a['R'] is None for a in arms) else sum(a['R'] for a in arms)/2
            value = {**classify(val, band, old['expected_class']),
                     'required_members': len(old['members']),
                     'missing_members': sorted({x for a in arms for x in a['missing_members']}),
                     'within_protein_range': None, 'equal_mean_of_complete_arm_medians': True}
        groups.append({**old, 'methods': {**old['methods'], METHOD: value}})
    triples = [{**old, 'methods': {**old['methods'], METHOD: strict_summary(
        old['members'], index, METHOD, band, old['expected_class'])}} for old in src['triples']]
    subsets = {'all225': rows, 'La100': [r for r in rows if r['source_conditioning_metal'] == 'La'],
               'Ca125': [r for r in rows if r['source_conditioning_metal'] == 'Ca'],
               **{p: [r for r in groups if r['pool'] == p] for p in ('La4', 'Ca5', 'balanced')},
               'La100_triples': triples}
    methods = ('context_composite', 'context_union', METHOD, ADAPTIVE)
    matched = {}
    for name, subset in subsets.items():
        common = [r for r in subset if all(r['methods'][m]['R'] is not None
                                         for m in (METHOD, ADAPTIVE))]
        matched[name] = {'common': len(common), 'counts': summarize_rows(common, (METHOD, ADAPTIVE)),
                         'static_to_adaptive': dict(Counter(r['methods'][METHOD]['outcome']+'->'+
                             r['methods'][ADAPTIVE]['outcome'] for r in common))}
    result = {'protocol_id': PROTOCOL, 'transfer': record(transfer), 'reference': record(reference),
              'rows': rows, 'pools': groups, 'triples': triples, 'bands': bands,
              'counts': {k: summarize_rows(v, methods) for k, v in subsets.items()},
              'matched': matched, 'new_molecular_calls': 0, 'production_changed': False,
              'interpretation': 'Consumed development sources; static and adaptive have separate canonical-only bands.',
              'implementation': record(__file__)}
    write_new(output, result)
    return result['counts']['all225']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='op', required=True)
    for op, fields in {'calibrate': ('qualification', 'output'),
                       'compare': ('transfer', 'reference', 'output')}.items():
        command = sub.add_parser(op)
        for field in fields:
            command.add_argument('--'+field, type=Path, required=True)
    args = vars(parser.parse_args())
    print(json.dumps(globals()[args.pop('op')](**args), indent=2))
