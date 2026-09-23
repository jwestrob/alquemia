"""Replay one declared reduced adaptive geometry pool; no molecular calls."""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from nikasha_pool import SCALED_PROTOCOL, choose_rows
from nikasha_pool_compare import (VARIANTS, load_results, original_reference,
    calibration_rows, extrema_reference, old_bands, method_record, aggregate)
from accommodation_fold_proposals import summarize_rows

PROTOCOL = 'Nikasha_origin_and_two_adaptive_candidates_replay_v1'
CANDIDATES = ('origin', 'adaptive_Ca', 'adaptive_La')
REFERENCE_ID = 'Nikasha_minimal_adaptive_canonical25_reference_v1'


def reduced_case(case):
    matrix = {'Ca': {}, 'La': {}}
    for candidate in CANDIDATES:
        representative = case.get('aliases', {}).get(candidate, {}).get('representative', candidate)
        for metal in matrix:
            cell = case.get('matrix', {}).get(metal, {}).get(representative)
            if cell is not None:
                matrix[metal][candidate] = cell
    pool = choose_rows(matrix, CANDIDATES)
    if case['status'] != 'prepared' and pool['status'] == 'available':
        raise InvalidArtifact('unprepared source unexpectedly has a complete matrix')
    if pool['status'] == 'available':
        # Removing candidates cannot lower a mathematical row minimum.
        for metal in matrix:
            old = case['pool']['rows'][metal]
            new = pool['rows'][metal]
            a = old['work_from_origin_kcal_mol'][old['mathematical_candidate']]['composite_kcal_mol']
            b = new['work_from_origin_kcal_mol'][new['mathematical_candidate']]['composite_kcal_mol']
            if b < a - 1e-8:
                raise InvalidArtifact('subset lowered a mathematical endpoint minimum')
    return {**case, 'full_pool': case['pool'], 'pool': pool,
            'declared_candidates': list(CANDIDATES)}


def bundle(path, population):
    b = load_results([path])
    if b['protocol_id'] != SCALED_PROTOCOL or len(b['manifests']) != 1:
        raise InvalidArtifact('one completed scaled-angular archive required')
    if b['manifests'][0]['population'] != population:
        raise InvalidArtifact('wrong declared source population')
    b['rows'] = {cid: reduced_case(row) for cid, row in b['rows'].items()}
    return b


def calibrate(collection, agreement, output):
    b = bundle(collection, 'scaled30')
    pin, original = original_reference(b)
    if len(b['rows']) != 30:
        raise InvalidArtifact('original30 denominator differs')
    variants = {v: extrema_reference(calibration_rows(b, original, v), v, REFERENCE_ID)
                for v in VARIANTS}
    result = {'protocol_id': PROTOCOL, 'candidate_ids': list(CANDIDATES),
        'reference_id': REFERENCE_ID, 'variants': variants, 'signature': b['signature'],
        'collection': record(collection), 'original_reference': pin,
        'old_frozen_bands': old_bands(original), 'agreement': record(agreement),
        'implementation': record(__file__), 'frozen_at_UTC': datetime.now(timezone.utc).isoformat(),
        'numerical_provenance': b['numerical_provenance'],
        'noncanonical_folds_used_for_calibration': False, 'new_molecular_calls': 0,
        'rule': 'original25 class extrema; unchanged minimum gap; no crystal or PLM calibration'}
    write_new(output, result)
    return {v: {k: variants[v][k] for k in ('status', 'bands', 'gap_model_kcal_mol')}
            for v in VARIANTS}


def checked_reference(path, b):
    r = read_json(path)
    if (r['protocol_id'] != PROTOCOL or r['candidate_ids'] != list(CANDIDATES)
            or r['signature'] != b['signature'] or r['noncanonical_folds_used_for_calibration']):
        raise InvalidArtifact('reference method/population policy differs')
    verify(r['agreement']); verify(r['collection']); verify(r['implementation'])
    original = read_json(verify(r['original_reference']))
    designated = {x['case_id']: x for x in original['rows']
                  if x['representation'] == 'context' and x['role'] == 'calibration'}
    for v in VARIANTS:
        rows = r['variants'][v]['rows']
        if (len(rows) != 25 or {x['case_id'] for x in rows} != set(designated)
                or any(x['expected_class'] != designated[x['case_id']]['expected_class'] for x in rows)
                or extrema_reference(rows, v, REFERENCE_ID) != r['variants'][v]):
            raise InvalidArtifact('canonical-only frozen reference differs')
    return r


def compare(collection, reference, full_comparison, output):
    b = bundle(collection, 'scaled225'); ref = checked_reference(reference, b)
    prior = read_json(full_comparison)
    if prior['collections'] != [record(collection)] or len(b['rows']) != 225:
        raise InvalidArtifact('full comparison source differs')
    if {r['case_id'] for r in prior['rows']} != set(b['rows']):
        raise InvalidArtifact('all225 source identities required')
    bands = {'released': prior['bands']['context_composite'],
             'full_adaptive': prior['bands']['pool_operational_new'],
             'minimal_old': ref['old_frozen_bands'],
             'minimal_full_bands': prior['bands']['pool_operational_new'],
             **{'minimal_' + v: ref['variants'][v]['bands'] for v in VARIANTS}}
    rows = []
    for old in prior['rows']:
        case = b['rows'][old['case_id']]; pool = case['pool']
        methods = {'released': old['methods']['context_composite'],
                   'full_adaptive': old['methods']['pool_operational_new']}
        for method in ('minimal_old', 'minimal_full_bands', 'minimal_mathematical', 'minimal_operational'):
            v = 'mathematical' if method == 'minimal_mathematical' else 'operational'
            value = pool[v]['composite_R_model_kcal_mol'] if pool[v] else None
            methods[method] = method_record(value, bands[method], old['expected_class'])
        rows.append({**old, 'methods': methods, 'minimal_pool': pool,
                     'full_pool': case['full_pool'], 'unavailable_reason': case.get('reason'),
                     'parent_status': case['status']})
    index = {r['case_id']: r for r in rows}
    pools, triples = [], []
    for group in prior['pools']:
        if group['pool'] == 'balanced':
            arm = {z: next(x for x in pools if x['root_case_id'] == group['root_case_id']
                          and x['pool'] == z) for z in ('La4', 'Ca5')}
            methods = {}
            for m, limits in bands.items():
                a, c = [arm[z]['methods'][m]['R'] for z in ('La4', 'Ca5')]
                value = (a + c)/2 if a is not None and c is not None else None
                methods[m] = {**method_record(value, limits, group['expected_class']),
                              'Ca5_minus_La4': c - a if value is not None else None}
        else:
            methods = {m: aggregate(group['members'], index, m, limits, group['expected_class'])
                       for m, limits in bands.items()}
        pools.append({**group, 'methods': methods})
    for group in prior['triples']:
        triples.append({**group, 'methods': {m: aggregate(group['members'], index, m, limits,
                        group['expected_class']) for m, limits in bands.items()}})
    subsets = {'all225': rows, 'La100': [r for r in rows if r['source_conditioning_metal'] == 'La'],
               'Ca125': [r for r in rows if r['source_conditioning_metal'] == 'Ca'],
               'all100_triples': triples,
               **{p: [r for r in pools if r['pool'] == p] for p in ('La4', 'Ca5', 'balanced')}}
    counts = {name: summarize_rows(values, bands) for name, values in subsets.items()}
    matched = {}
    for name, values in subsets.items():
        matched[name] = {}
        for comparator in ('released', 'full_adaptive'):
            common = [r for r in values if all(r['methods'][m]['decision'] != 'unavailable'
                      for m in (comparator, 'minimal_operational'))]
            matched[name][comparator] = {'common': len(common),
                'counts': summarize_rows(common, (comparator, 'minimal_operational')),
                'transitions': dict(Counter(r['methods'][comparator]['outcome'] + '->' +
                    r['methods']['minimal_operational']['outcome'] for r in common))}
    result = {'protocol_id': PROTOCOL, 'candidate_ids': list(CANDIDATES),
        'collection': record(collection), 'reference': record(reference),
        'full_comparison': record(full_comparison), 'implementation': record(__file__),
        'bands': bands, 'rows': rows, 'pools': pools, 'triples': triples,
        'counts': counts, 'matched': matched, 'new_molecular_calls': 0,
        'production_changed': False, 'new_numerical_qualification': False,
        'nominal_distinct_candidate_GFNs': {'minimal': 12, 'full': 20},
        'interpretation': 'consumed structural development data; call counts are not measured latency'}
    write_new(output, result)
    return {'all225': counts['all225'], 'matched225': matched['all225']}


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='op', required=True)
    for op in ('calibrate', 'compare'):
        a = sub.add_parser(op)
        for name in ('collection', 'output', *({'calibrate': ('agreement',),
                    'compare': ('reference', 'full-comparison')}[op])):
            a.add_argument('--' + name, required=True, type=Path)
    args = vars(p.parse_args()); op = args.pop('op')
    print(json.dumps(globals()[op](**args), indent=2))


if __name__ == '__main__':
    main()
