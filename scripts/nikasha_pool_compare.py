"""Canonical-only references and fixed comparisons for actual shared-pool results."""
from __future__ import annotations
import argparse
from collections import Counter
from itertools import combinations
from datetime import datetime, timezone
import json
from pathlib import Path
import statistics

from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from accommodation_fold_proposals import outcome, summarize_rows, strict_summary
from accommodation_folds_compare import source_groups, decision
from compact_solvation_compare import MINIMUM_CALIBRATION_GAP
from nikasha_pool import PROTOCOL, SETTINGS, ADAPTIVE_PROTOCOL, ADAPTIVE_SETTINGS, SCALED_PROTOCOL, choose_rows

VARIANTS = ('mathematical', 'operational')
REFERENCE_ID = 'Nikasha_common_geometry_canonical25_reference_v1'
ADAPTIVE_REFERENCE_ID = 'Nikasha_adaptive_angular_common_geometry_canonical25_reference_v1'
SCALED_REFERENCE_ID = 'Nikasha_scaled_angular_common_geometry_canonical25_reference_v1'
PROTOCOLS = {
    PROTOCOL: (SETTINGS, REFERENCE_ID, {'pilot4', 'remaining26'}),
    ADAPTIVE_PROTOCOL: (ADAPTIVE_SETTINGS, ADAPTIVE_REFERENCE_ID, {'adaptive4', 'adaptive26'}),
    SCALED_PROTOCOL: (ADAPTIVE_SETTINGS, SCALED_REFERENCE_ID, {'scaled30'}),
}


def numerical_provenance(data, manifest, collection):
    """Retain qualified iteration-ceiling differences without treating them as a new Hamiltonian."""
    ceiling = manifest.get('GFN2_maxiter', 125)
    if ceiling not in (125, 500): raise InvalidArtifact('unsupported adaptive numerical ceiling')
    policy = manifest.get('numerical_policy_id', 'native_GFN2_MaxIter125_unchanged_convergence_v1')
    if policy != f'native_GFN2_MaxIter{ceiling}_unchanged_convergence_v1':
        raise InvalidArtifact('adaptive numerical policy differs')
    recovered = data.get('recovery_overlay_applied', False)
    qualification_pin = data.get('numerical_diagnostic') if recovered else manifest.get('numerical_qualification')
    if ceiling == 500 or recovered:
        if qualification_pin is None: raise InvalidArtifact('iteration ceiling lacks actual qualification')
        qualification = read_json(verify(qualification_pin)); qm = read_json(verify(qualification['manifest']))
        if (qualification['protocol_id'] != 'adaptive_pool_native_GFN2_MaxIter500_diagnostic_v1' or
                qualification['complete'] != 3 or qualification['denominator'] != 3 or
                not qualification['both_controls_pass'] or not qualification['formerly_failed_cell_complete'] or
                qm['orca'] != manifest['orca'] or qm['control_tolerance_kcal_mol'] != .05 or
                any(r['status'] != 'complete' or not r['parameters_identical'] or
                    r['charge_sanity_status'] != 'pass' for r in qualification['rows'])):
            raise InvalidArtifact('incompatible adaptive numerical qualification')
    if recovered:
        primary = read_json(verify(data['recovered_from']))
        primary_rows = {r['case_id']: r for r in primary['cases']}
        if (data.get('numerical_policy_id') != 'primary_native_GFN2_with_explicit_MaxIter500_completion_v1' or
                primary['protocol_id'] != ADAPTIVE_PROTOCOL or primary['manifest'] != data['manifest'] or
                data['additional_GFN2_attempts'] != 3 or
                data['primary_GFN2_complete'] != primary['GFN2_complete'] or
                set(primary_rows) != {r['case_id'] for r in data['cases']} or
                any(r.get('primary_pool') != primary_rows[r['case_id']]['pool'] for r in data['cases'])):
            raise InvalidArtifact('adaptive recovery lost its primary failed result')
    return {'collection': record(collection), 'manifest_GFN2_maxiter': ceiling,
            'manifest_numerical_policy_id': policy,
            'collection_numerical_policy_id': data.get('numerical_policy_id', policy),
            'qualification': qualification_pin, 'recovery_overlay_applied': recovered,
            'recovered_from': data.get('recovered_from'),
            'primary_GFN2_complete': data.get('primary_GFN2_complete', data.get('GFN2_complete')),
            'additional_GFN2_attempts': data.get('additional_GFN2_attempts', 0)}


def load_results(collections):
    """Check saved matrix algebra; do not rerun the source output parsers or models."""
    rows, manifests, pins, signature, protocol, numerical = {}, [], [], None, None, []
    for path in collections:
        data = read_json(path); m = read_json(verify(data['manifest']))
        current = data['protocol_id']
        if (current not in PROTOCOLS or m['protocol_id'] != current or
                m['settings'] != PROTOCOLS[current][0]):
            raise InvalidArtifact('shared-pool protocol differs')
        if protocol is not None and protocol != current: raise InvalidArtifact('mixed pool protocols')
        protocol = current
        if protocol in (ADAPTIVE_PROTOCOL, SCALED_PROTOCOL):
            numerical.append(numerical_provenance(data, m, path))
        sig = {k: m[k] for k in ('model', 'software', 'orca', 'settings')}
        if protocol == SCALED_PROTOCOL:
            sig.update({k: m[k] for k in ('proposal_protocol_id', 'proposal_settings')})
            proposer = read_json(verify(read_json(verify(m['adaptive_proposals']))['manifest']))
            if m['proposal_protocol_id'] != proposer['protocol_id'] or m['proposal_settings'] != proposer['settings']:
                raise InvalidArtifact('scaled proposal source policy changed')
        if signature is not None and sig != signature: raise InvalidArtifact('mixed pool methods or software')
        signature = sig
        if len(data['cases']) != data['denominator'] or {c['case_id'] for c in data['cases']} != set(m['declared_case_ids']):
            raise InvalidArtifact('pool population differs from its manifest')
        source = read_json(verify(m['source'])); originals = {c['case_id']: c for c in source['cases']}
        for c in data['cases']:
            cid = c['case_id']
            if cid in rows: raise InvalidArtifact('duplicate pool identity across collections: ' + cid)
            if c['old_result'] != originals[cid]: raise InvalidArtifact('archived own-proposal result changed')
            if c['status'] == 'prepared':
                calculated = choose_rows(c['matrix'], [q['id'] for q in c['candidates']])
                if calculated != c['pool']: raise InvalidArtifact('saved matrix selection/score algebra differs')
            elif c['pool']['status'] != 'unavailable' or any(c['pool'].get(v) is not None for v in VARIANTS):
                raise InvalidArtifact('unsupported source incorrectly scored')
            rows[cid] = {**c, 'collection': record(path)}
        manifests.append(m); pins.append(record(path))
    if not manifests: raise InvalidArtifact('no actual pool collections provided')
    return {'rows': rows, 'manifests': manifests, 'collections': pins, 'signature': signature,
            'protocol_id': protocol, 'numerical_provenance': numerical}


def original_reference(bundle):
    references = []
    for m in bundle['manifests']:
        parent = read_json(verify(m['source_manifest']))
        references.append(parent['frozen_comparison'])
    if any(p != references[0] for p in references): raise InvalidArtifact('original reference differs across sources')
    return references[0], read_json(verify(references[0]))


def old_bands(original):
    b = original['calibration']['context']['bands']
    return {'Ca_max': b['Ca_supported_max_R_model_kcal_mol'], 'La_min': b['La_supported_min_R_model_kcal_mol']}


def calibration_rows(bundle, original, variant):
    designated = [r for r in original['rows'] if r['representation'] == 'context' and r['role'] == 'calibration']
    if len(designated) != 25 or len({r['case_id'] for r in designated}) != 25:
        raise InvalidArtifact('released canonical25 membership differs')
    result = []
    for src in designated:
        case = bundle['rows'].get(src['case_id']); score = None
        if case:
            old = case['old_result']
            if (old['expected_class'], old['label_scope']) != (src['expected_class'], src['label_scope']):
                raise InvalidArtifact('canonical biological label/stratum changed')
            score = case['pool'][variant]
        result.append({'case_id': src['case_id'], 'expected_class': src['expected_class'],
                       'biological_group': src['biological_group'], 'label_scope': src['label_scope'],
                       'collection': case['collection'] if case else None,
                       'status': 'available' if score is not None else 'unavailable',
                       'R_model_kcal_mol': score['composite_R_model_kcal_mol'] if score else None,
                       'selected_candidates': {z: case['pool']['rows'][z][variant + '_candidate'] for z in ('Ca', 'La')} if score else None})
    return result


def extrema_reference(rows, variant, reference_id=REFERENCE_ID):
    available = [r for r in rows if r['R_model_kcal_mol'] is not None]
    band = gap = extrema = None; status = 'unavailable_calibration_member'
    if len(rows) != 25: raise InvalidArtifact('calibration requires all25 designated members')
    if len(available) == 25:
        values = {z: [r['R_model_kcal_mol'] for r in rows if r['expected_class'] == z] for z in ('Ca', 'La')}
        if not all(values.values()): raise InvalidArtifact('missing canonical class')
        extrema = {'Ca_max': max(values['Ca']), 'La_min': min(values['La'])}
        gap = extrema['La_min'] - extrema['Ca_max']
        status = 'available' if gap > MINIMUM_CALIBRATION_GAP else 'unsupported_class_separation'
        if status == 'available': band = extrema
    return {'reference_id': reference_id + '_' + variant, 'variant': variant, 'status': status,
            'calibration_denominator': 25, 'available_calibration': len(available), 'bands': band,
            'class_extrema': extrema, 'gap_model_kcal_mol': gap,
            'minimum_gap_model_kcal_mol': MINIMUM_CALIBRATION_GAP, 'rows': rows}


def calibrate(collections, agreement, output):
    bundle = load_results(collections); original_pin, original = original_reference(bundle)
    protocol = bundle['protocol_id']; _, reference_id, populations = PROTOCOLS[protocol]
    if len(bundle['manifests']) != len(populations) or {m['population'] for m in bundle['manifests']} != populations:
        raise InvalidArtifact('reference requires the declared ' + ' + '.join(sorted(populations)) + ' original sources')
    parents = [read_json(verify(m['source_manifest'])) for m in bundle['manifests']]
    expected = {c['case_id'] for c in parents[0]['cases']}
    if len(expected) != 30 or set(bundle['rows']) != expected or any({c['case_id'] for c in p['cases']} != expected for p in parents):
        raise InvalidArtifact('original30 source union differs; folds may not enter calibration')
    variants = {v: extrema_reference(calibration_rows(bundle, original, v), v, reference_id) for v in VARIANTS}
    result = {'protocol_id': protocol, 'reference_id': reference_id, **bundle['signature'],
              'variants': variants, 'collections': bundle['collections'], 'original_reference': original_pin,
              'old_frozen_bands': old_bands(original), 'agreement': record(agreement), 'implementation': record(__file__),
              'frozen_at_UTC': datetime.now(timezone.utc).isoformat(),
              'rule': 'released canonical25 class extrema and unchanged minimum gap; each score variant separately',
              'noncanonical_folds_used_for_calibration': False, 'crystals_or_PLM_used_for_calibration': False,
              'new_molecular_calls': 0, 'production_changed': False,
              'interpretation': 'developmental calibration, not independent validation'}
    if protocol in (ADAPTIVE_PROTOCOL, SCALED_PROTOCOL): result['numerical_provenance'] = bundle['numerical_provenance']
    write_new(output, result)
    return {v: {k: variants[v][k] for k in ('status', 'bands', 'gap_model_kcal_mol', 'available_calibration')} for v in VARIANTS}


def checked_reference(path, bundle):
    reference = read_json(path)
    protocol = bundle['protocol_id']; reference_id = PROTOCOLS[protocol][1]
    if (reference['protocol_id'] != protocol or reference['reference_id'] != reference_id or
            reference['noncanonical_folds_used_for_calibration'] or reference['crystals_or_PLM_used_for_calibration']):
        raise InvalidArtifact('incompatible or contaminated pool reference')
    if any(reference[k] != bundle['signature'][k] for k in bundle['signature']):
        raise InvalidArtifact('pool reference method differs')
    original = read_json(verify(reference['original_reference']))
    designated = {r['case_id']: r for r in original['rows'] if r['representation'] == 'context' and r['role'] == 'calibration'}
    if len(designated) != 25: raise InvalidArtifact('reference canonical membership differs')
    for v in VARIANTS:
        saved = reference['variants'][v]
        if len(saved['rows']) != 25 or {r['case_id'] for r in saved['rows']} != set(designated):
            raise InvalidArtifact('saved reference includes different calibration members')
        if any(r['expected_class'] != designated[r['case_id']]['expected_class'] for r in saved['rows']):
            raise InvalidArtifact('saved canonical labels differ')
        # Integrity check of the frozen artifact, not a new fit on transfer data.
        if extrema_reference(saved['rows'], v, reference_id) != saved:
            raise InvalidArtifact('saved reference extrema/bands differ from its canonical rows')
    numerical = []
    for pin in reference['collections']:
        source = verify(pin)
        if protocol in (ADAPTIVE_PROTOCOL, SCALED_PROTOCOL):
            data = read_json(source)
            numerical.append(numerical_provenance(data, read_json(verify(data['manifest'])), source))
    if protocol in (ADAPTIVE_PROTOCOL, SCALED_PROTOCOL) and reference.get('numerical_provenance') != numerical:
        raise InvalidArtifact('saved adaptive numerical provenance differs')
    return reference


def call(value, bands):
    return decision(value, bands) if bands is not None else 'unavailable'


def method_record(value, bands, expected):
    label = call(value, bands)
    return {'R': value, 'decision': label, 'outcome': outcome(label, expected),
            'score_status': 'available' if value is not None else 'unavailable',
            'reference_status': 'available' if bands is not None else 'unavailable'}


def inspect(collections, output, reference=None):
    bundle = load_results(collections); original_pin, original = original_reference(bundle)
    ref = checked_reference(reference, bundle) if reference else None; bands = old_bands(original); rows = []
    for cid, case in bundle['rows'].items():
        old = case['old_result']; pool = case['pool']
        static = old['R0']['composite_R_model_kcal_mol'] if old['R0'] else None
        own = old['R_selected']['composite_R_model_kcal_mol'] if old['R_selected'] else None
        r = {'case_id': cid, 'expected_class': old['expected_class'], 'pool_status': pool['status'],
             'static_R': static, 'own_proposal_R': own, 'variants': {}, 'collection': case['collection']}
        if bundle['protocol_id'] in (ADAPTIVE_PROTOCOL, SCALED_PROTOCOL):
            r.update(primary_pool=case.get('primary_pool', pool), prior_pool=case['prior_pool'])
        for v in VARIANTS:
            value = pool[v]['composite_R_model_kcal_mol'] if pool[v] else None
            new = ref['variants'][v]['bands'] if ref else None
            r['variants'][v] = {'R': value, 'old_band_decision': call(value, bands), 'new_band_decision': call(value, new),
                'delta_from_static': value - static if value is not None and static is not None else None,
                'delta_from_own_proposal': value - own if value is not None and own is not None else None,
                'selected_candidates': {z: pool['rows'][z][v + '_candidate'] for z in ('Ca', 'La')} if pool['rows'] else None}
        rows.append(r)
    result = {'protocol_id': bundle['protocol_id'], 'collections': bundle['collections'], 'reference': record(reference) if reference else None,
              'original_reference': original_pin, 'rows': rows, 'new_molecular_calls': 0,
              'interpretation': 'raw development comparison; unknown PLM labels remain unknown'}
    if bundle['protocol_id'] == ADAPTIVE_PROTOCOL: result['numerical_provenance'] = bundle['numerical_provenance']
    write_new(output, result)
    lines = ['# Actual shared-pool contrast changes', '', '| Case | Variant | R | Δ from own proposal | Old-band transfer | New-band call |', '|---|---|---:|---:|---|---|']
    for row in rows:
        for v, d in row['variants'].items():
            values = ['unavailable' if d[k] is None else f'{d[k]:.6f}' for k in ('R', 'delta_from_own_proposal')]
            lines.append('| ' + row['case_id'] + ' | ' + v + ' | ' + ' | '.join(values) + ' | ' + d['old_band_decision'] + ' | ' + d['new_band_decision'] + ' |')
    lines += ['', 'Old-band calls are transfer checks. New calls are unavailable until a compatible variant-specific reference exists. No affinity interpretation, PLM truth labels or population weights are inferred.', '']
    with Path(output).with_suffix('.md').open('x') as f: f.write('\n'.join(lines))
    return {'rows': len(rows), 'pool_available': sum(r['pool_status'] == 'available' for r in rows)}


def aggregate(members, index, method, bands, expected):
    if bands is not None: return strict_summary(members, index, method, bands, expected)
    values = [index[c]['methods'][method]['R'] for c in members]
    complete = all(v is not None for v in values)
    return {**method_record(statistics.median(values) if complete else None, None, expected),
            'required_members': len(members), 'missing_members': [c for c, v in zip(members, values) if v is None],
            'within_protein_range': max(values) - min(values) if complete else None}


def compare(collections, reference, prior_comparison, output):
    bundle = load_results(collections); ref = checked_reference(reference, bundle); prior = read_json(prior_comparison)
    population = {PROTOCOL: 'primary225', SCALED_PROTOCOL: 'scaled225'}.get(bundle['protocol_id'])
    if population is None or len(bundle['manifests']) != 1 or bundle['manifests'][0]['population'] != population:
        raise InvalidArtifact('fold comparison requires the full primary225 pool collection')
    pm = read_json(verify(prior['manifest'])); source = read_json(verify(pm['sources'])); groups = source_groups(source)
    if bundle['manifests'][0]['source'] != prior['selection']:
        raise InvalidArtifact('pool and prior comparison describe different original proposal results')
    expected_ids = {c['case_id'] for g in groups.values() for c in g if c['primary_evaluation_pool']}
    if set(bundle['rows']) != expected_ids: raise InvalidArtifact('all225 primary identities must remain')
    bands = dict(prior['bands'])
    for v in VARIANTS:
        bands['pool_' + v + '_old'] = ref['old_frozen_bands']
        bands['pool_' + v + '_new'] = ref['variants'][v]['bands']
    rows = []
    for old in prior['rows']:
        c = bundle['rows'][old['case_id']]; row = {**old, 'methods': dict(old['methods'])}
        if c['old_result']['expected_class'] != old['expected_class']: raise InvalidArtifact('fold label changed')
        for v in VARIANTS:
            value = c['pool'][v]['composite_R_model_kcal_mol'] if c['pool'][v] else None
            for band_kind in ('old', 'new'):
                key = 'pool_' + v + '_' + band_kind
                row['methods'][key] = method_record(value, bands[key], row['expected_class'])
        row['pool_status'] = c['pool']['status']; rows.append(row)
    index = {r['case_id']: r for r in rows}; pools, triples = [], []
    for parent, members in groups.items():
        expected = members[0]['expected_class']
        arms = {z: [r['case_id'] for r in members if r['primary_evaluation_pool'] and r['source_conditioning_metal'] == z] for z in ('La', 'Ca')}
        summaries = {z: {m: aggregate(arms[z], index, m, b, expected) for m, b in bands.items()} for z in arms}
        for z, name in (('La', 'La4'), ('Ca', 'Ca5')):
            pools.append({'root_case_id': parent, 'expected_class': expected, 'pool': name, 'members': arms[z], 'methods': summaries[z]})
        balanced = {}
        for method, b in bands.items():
            a, c = summaries['La'][method]['R'], summaries['Ca'][method]['R']
            value = (a + c) / 2 if a is not None and c is not None else None
            balanced[method] = {**method_record(value, b, expected), 'Ca5_minus_La4': c - a if value is not None else None}
        pools.append({'root_case_id': parent, 'expected_class': expected, 'pool': 'balanced', 'members': arms['La'] + arms['Ca'], 'methods': balanced})
        for ids in combinations(arms['La'], 3):
            triples.append({'root_case_id': parent, 'expected_class': expected, 'members': list(ids),
                            'methods': {m: aggregate(ids, index, m, b, expected) for m, b in bands.items()}})
    if (len(rows), len(pools), len(triples)) != (225, 75, 100): raise InvalidArtifact('comparison denominators differ')
    subsets = {'all225': rows, 'La100': [r for r in rows if r['source_conditioning_metal'] == 'La'],
               'Ca125': [r for r in rows if r['source_conditioning_metal'] == 'Ca'], 'all100_triples': triples,
               **{p: [r for r in pools if r['pool'] == p] for p in ('La4', 'Ca5', 'balanced')}}
    counts = {name: summarize_rows(subset, bands) for name, subset in subsets.items()}; matched = {}
    for name, subset in subsets.items():
        matched[name] = {}
        for pool_method in [m for m in bands if m.startswith('pool_')]:
            matched[name][pool_method] = {}
            for prior_method in prior['bands']:
                common = [r for r in subset if all(r['methods'][m]['R'] is not None and r['methods'][m]['decision'] != 'unavailable' for m in (pool_method, prior_method))]
                matched[name][pool_method][prior_method] = {'declared': len(subset), 'common': len(common),
                    'counts': summarize_rows(common, (pool_method, prior_method)),
                    'outcome_transitions': dict(Counter(r['methods'][prior_method]['outcome'] + '->' + r['methods'][pool_method]['outcome'] for r in common))}
    result = {'protocol_id': bundle['protocol_id'], 'collections': bundle['collections'], 'reference': record(reference),
              'prior_comparison': record(prior_comparison), 'sources': pm['sources'], 'bands': bands,
              'rows': rows, 'pools': pools, 'triples': triples, 'counts': counts, 'matched': matched,
              'implementation': record(__file__), 'new_molecular_calls': 0, 'thresholds_fitted_on_folds': False,
              'production_changed': False, 'interpretation': 'consumed protein structural repeats; all pools/triples require every member'}
    write_new(output, result)
    lines = ['# Shared-pool fixed transfer comparison', '',
             'All cases retained. Mathematical and operational variants use separate canonical-only references; old-band transfer remains visible. Structural repeats are not independent biological labels.', '',
             '| Population | Method | Correct | Wrong | Inconclusive | Unavailable | Total |', '|---|---|---:|---:|---:|---:|---:|']
    for population, methods in counts.items():
        for method, tally in methods.items():
            lines.append('| ' + population + ' | ' + method + ' | ' + ' | '.join(str(tally[k]) for k in ('correct', 'wrong', 'inconclusive', 'unavailable', 'denominator')) + ' |')
    lines += ['', 'JSON retains raw contrasts, all strict members, within-protein ranges, both pool choices, and matched comparison denominators. A missing reference makes the decision unavailable without discarding a computed raw contrast.', '']
    with Path(output).with_suffix('.md').open('x') as f: f.write('\n'.join(lines))
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='operation', required=True)
    for name in ('inspect', 'calibrate', 'compare'):
        p = sub.add_parser(name); p.add_argument('--collections', nargs='+', required=True); p.add_argument('--output', required=True)
        if name == 'calibrate': p.add_argument('--agreement', required=True)
        if name == 'inspect': p.add_argument('--reference')
        if name == 'compare':
            p.add_argument('--reference', required=True); p.add_argument('--prior-comparison', required=True)
    args = vars(parser.parse_args()); op = args.pop('operation'); print(json.dumps(globals()[op](**args)))


if __name__ == '__main__': main()
