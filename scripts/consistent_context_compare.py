"""Actual component collection and canonical-only reference for union contexts."""
from __future__ import annotations
import argparse
from collections import Counter
import datetime
import json
from pathlib import Path

from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from consistent_context import PROTOCOL, PILOT_GROUPS, reusable_state, validate_solvent
from compact_solvation_compare import native_endpoint, mix_pair, MINIMUM_CALIBRATION_GAP
import compact_solvation as solvent
from mace_hybrid import EV_TO_KCAL
from accommodation_folds_compare import decision, source_groups, pool_summary, COMPONENTS
from accommodation_folds_triples import outcome


def collect(stage, output):
    stage = Path(stage).resolve(); ready = read_json(stage/'READY.json')
    inv = read_json(verify(ready['inventory'])); mp = verify(ready['MACE_manifest']); mm = read_json(mp)
    sp = verify(ready['GFN2_manifest']); sm = read_json(sp); validate_solvent(sp)
    if sm['inventory'] != ready['inventory'] or mm['preparation'] != ready['preparation']:
        raise InvalidArtifact('stage source manifests differ')
    rows = []
    for c in inv['cases']:
        row = {k: v for k, v in c.items() if k not in ('representations', 'old_context')}
        native = {}; low = {}; new_native = 0; new_low = 0
        for z in ('Ca', 'La'):
            tid = c['case_id']+'__context__'+z
            try:
                rp = verify(mm['reused'][tid]['receipt']) if tid in mm['reused'] else mp.parent/'results'/tid/'result.json'
                e = native_endpoint(read_json(rp), mm['model'])
                if not reusable_state(c['representations']['context']['endpoints'][z], e):
                    raise InvalidArtifact('actual MACE result has a different source state')
                native[z] = {**e, 'status': 'complete', 'native_parser_status': e['status'], 'reused': tid in mm['reused']}
                new_native += int(tid not in mm['reused'])
            except (OSError, InvalidArtifact, KeyError) as exc:
                native[z] = {'status': 'unavailable', 'reason': str(exc), 'native_MACE_energy_eV': None}
            low[z] = {}
            for medium in ('vacuum', 'alpb'):
                task = next(t for t in sm['all_tasks'] if (t['case_id'], t['metal'], t['medium']) == (c['case_id'], z, medium))
                reused = task['task_id'] in sm['reused']
                pin = sm['reused'].get(task['task_id']) or solvent.completed(sp, task['task_id'])
                if pin:
                    try:
                        # Audit the original actual execution, not a relabeled copy.
                        original = read_json(verify(pin['manifest']))
                        original_task = next(t for t in original['all_tasks'] if t['task_id'] == pin['task_id'])
                        parsed = solvent.diagnostics(pin, original_task)
                        low[z][medium] = {'status': 'complete', **pin, 'reused': reused,
                                          'current_task_id': task['task_id'], 'component_audit': parsed}
                        new_low += int(not reused)
                    except (OSError, InvalidArtifact, KeyError) as exc:
                        low[z][medium] = {'status': 'unavailable', 'reason': str(exc), 'energy_hartree': None, 'reused': reused}
                else:
                    op = Path(task['output_path']); rp = Path(str(op)+'.execution.json')
                    low[z][medium] = {'status': 'unavailable', 'energy_hartree': None, 'reused': False,
                                      'available_artifacts': [record(p) for p in (op, rp) if p.exists()]}
        native_ready = all(native[z]['status'] == 'complete' for z in native)
        complete = native_ready and all(e['status'] == 'complete' for r in low.values() for e in r.values())
        row.update(status='complete' if complete else 'unavailable', native_endpoints=native, solvent_endpoints=low,
                   native_R_model_kcal_mol=(native['Ca']['native_MACE_energy_eV']-native['La']['native_MACE_energy_eV'])*EV_TO_KCAL if native_ready else None,
                   composite_R_model_kcal_mol=None, solvation_delta_R_kcal_mol=None,
                   complete_fresh_MACE_calls=new_native, complete_fresh_GFN2_calls=new_low,
                   source_preparation=c['representations']['context']['preparation'],
                   changes_are_context_composition_cavity_and_preparation_not_pure_geometry=True)
        if complete:
            row.update(mix_pair({z: native[z]['native_MACE_energy_eV'] for z in native},
                       {z: low[z]['vacuum']['energy_hartree'] for z in low}, {z: low[z]['alpb']['energy_hartree'] for z in low}))
        rows.append(row)
    rows.extend({'case_id': f['case_id'], 'root_case_id': f['source']['root_case_id'],
                 'source_conditioning_metal': f['source']['source_conditioning_metal'],
                 'canonical_coordinate_match': f['source']['canonical_coordinate_match'],
                 'primary_evaluation_pool': f['source']['primary_evaluation_pool'],
                 'expected_class': f['source']['expected_class'], 'status': 'preparation_unavailable',
                 'reason': f['reason'], 'native_R_model_kcal_mol': None, 'composite_R_model_kcal_mol': None,
                 'solvation_delta_R_kcal_mol': None} for f in inv['failures'])
    result = {'protocol_id': PROTOCOL, 'stage': inv['stage'], 'ready': record(stage/'READY.json'),
              'inventory': ready['inventory'], 'MACE_manifest': ready['MACE_manifest'], 'GFN2_manifest': ready['GFN2_manifest'],
              'denominator': inv['denominator'], 'complete': sum(r['status'] == 'complete' for r in rows),
              'rows': rows, 'new_reference': None, 'model_or_default_changed': False,
              'analysis_implementation': {name: record(Path(__file__).with_name(name)) for name in
                   ('consistent_context_compare.py', 'consistent_context.py', 'compact_solvation.py', 'compact_solvation_compare.py')}}
    write_new(output, result)
    return {k: v for k, v in result.items() if k != 'rows'}


def reference(collection, output):
    result = read_json(collection)
    if result['protocol_id'] != PROTOCOL or result['stage'] != 'calibration28' or result['denominator'] != 28:
        raise InvalidArtifact('reference requires the declared canonical25/crystal3 collection')
    inv = read_json(verify(result['inventory'])); canonical = [c['case_id'] for c in inv['cases'] if c['role'] == 'calibration']
    # Include any explicitly unavailable canonical source in the strict denominator.
    canonical += [c['case_id'] for c in inv['failures'] if c['source']['canonical_coordinate_match']]
    if len(canonical) != 25 or len(set(canonical)) != 25: raise InvalidArtifact('original canonical25 membership changed')
    rows = [r for r in result['rows'] if r['case_id'] in canonical]
    ref = {'protocol_id': PROTOCOL, 'reference_id': 'Nikasha_consistent_context_canonical25_extrema_v1',
           'collection': record(collection), 'canonical_case_ids': sorted(canonical), 'denominator': 25,
           'available': sum(r['status'] == 'complete' for r in rows), 'status': 'unavailable_incomplete_canonical',
           'class_extrema': None, 'gap_model_kcal_mol': None, 'bands': None,
           'minimum_gap_model_kcal_mol': MINIMUM_CALIBRATION_GAP, 'transfers_used': False,
           'frozen_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'threshold_coefficients_fitted': False}
    if ref['available'] == 25:
        ca = [r['composite_R_model_kcal_mol'] for r in rows if r['expected_class'] == 'Ca']
        la = [r['composite_R_model_kcal_mol'] for r in rows if r['expected_class'] == 'La']
        if len(ca)+len(la) != 25 or not ca or not la: raise InvalidArtifact('unsupported canonical labels')
        gap = min(la)-max(ca)
        ref.update(class_extrema={'Ca_max': max(ca), 'La_min': min(la)}, gap_model_kcal_mol=gap,
                   status='available' if gap > MINIMUM_CALIBRATION_GAP else 'unavailable_overlapping_or_unresolved_classes',
                   class_spread={'Ca': max(ca)-min(ca), 'La': max(la)-min(la)},
                   ordered_crossclass_pairs=sum(a>b for a in la for b in ca), crossclass_pairs=len(ca)*len(la))
        if ref['status'] == 'available': ref['bands'] = {'Ca_max': max(ca), 'La_min': min(la)}
    write_new(output, ref)
    return ref


def report(collection, reference_path, static_reference, output):
    result = read_json(collection); ref = read_json(reference_path); static = read_json(static_reference)
    if ref['protocol_id'] != PROTOCOL: raise InvalidArtifact('reference protocol differs')
    b = static['calibration']['context']['bands']
    oldbands = {'Ca_max': b['Ca_supported_max_R_model_kcal_mol'], 'La_min': b['La_supported_min_R_model_kcal_mol']}
    table = []
    for row in result['rows']:
        value = row['composite_R_model_kcal_mol']
        table.append({k: row.get(k) for k in ('case_id', 'root_case_id', 'role', 'status', 'expected_class', 'native_R_model_kcal_mol',
                      'solvation_delta_R_kcal_mol', 'composite_R_model_kcal_mol')} |
                     {'old_band_transfer': decision(value, oldbands),
                      'own_reference_decision': decision(value, ref['bands']) if ref['bands'] else 'unavailable_reference'})
    summary = {}
    for role in sorted({r['role'] for r in table}, key=str):
        members = [r for r in table if r['role'] == role]; counts = Counter(r['own_reference_decision'] for r in members)
        summary[str(role)] = {'denominator': len(members), 'correct': sum(r['own_reference_decision'] == r['expected_class']+'-supported' for r in members),
                  'wrong': sum(r['own_reference_decision'] in ('La-supported', 'Ca-supported') and
                         r['own_reference_decision'] != r['expected_class']+'-supported' for r in members), 'decisions': dict(counts)}
    value = {'protocol_id': PROTOCOL, 'collection': record(collection), 'reference': record(reference_path),
             'static_reference': record(static_reference), 'summary': summary, 'rows': table}
    write_new(output, value)
    return {k: v for k, v in value.items() if k != 'rows'}


def transfer_compare(collection, reference_path, baseline, sources, output, scope='pilot36'):
    """Reuse frozen source groups and strict pool helper; no fitted transfer rules."""
    result = read_json(collection); ref = read_json(reference_path); old = read_json(baseline)
    src = read_json(sources); groups = source_groups(src)
    if ref['protocol_id'] != PROTOCOL or ref['status'] != 'available' or result['protocol_id'] != PROTOCOL:
        raise InvalidArtifact('valid same-protocol canonical reference required')
    inv = read_json(verify(result['inventory']))
    if inv['reference'] != record(reference_path):
        raise InvalidArtifact('transfer collection was prepared against another reference')
    roots = set(PILOT_GROUPS) if scope == 'pilot36' else set(groups) if scope == 'primary225' else set()
    if not roots: raise InvalidArtifact('undeclared transfer scope')
    selected = [r for r in src['cases'] if r['root_case_id'] in roots and r['primary_evaluation_pool']]
    if len(selected) != (36 if scope == 'pilot36' else 225): raise InvalidArtifact('declared source denominator differs')
    idx = {r['case_id']: r for r in result['rows']}; previous = {r['case_id']: r for r in old['rows']}
    if set(idx) != {r['case_id'] for r in selected}: raise InvalidArtifact('collection does not retain all declared transfer rows')
    rows = []
    for source in selected:
        row = idx[source['case_id']]; before = previous[source['case_id']]
        if row['expected_class'] != source['expected_class'] or before['expected_class'] != source['expected_class']:
            raise InvalidArtifact('original label changed')
        value = row['composite_R_model_kcal_mol']; call = decision(value, ref['bands'])
        methods = {k: dict(v) for k, v in before['methods'].items()}
        methods['context_union'] = {'R': value, 'decision': call, 'outcome': outcome(call, source['expected_class']),
                                   'old_band_transfer': decision(value, old['bands']['context_composite'])}
        original = methods['context_composite']['R']
        rows.append({k: source[k] for k in ('case_id', 'root_case_id', 'source_conditioning_metal', 'expected_class', 'biological_group')} |
                    {'methods': methods, 'union_status': row['status'],
                     'union_minus_static_R': None if value is None or original is None else value-original})
    # Existing strict helper requires every member and keeps all raw components.
    component_idx = {cid: {'case_id': cid, **{k: row.get(k) for k in COMPONENTS}} for cid, row in idx.items()}
    pool_bands = {'native': old['bands']['context_native'], 'composite': ref['bands']}
    pools = []; triples = []
    for parent in old['pools']:
        if parent['root_case_id'] not in roots: continue
        members = [next(r for r in selected if r['case_id'] == cid) for cid in parent['members']]
        if parent['pool'] in ('La4', 'Ca5'):
            summary = pool_summary(members, component_idx, pool_bands)
            summary['native_old_band_transfer'] = summary.pop('native_decision')
            summary['native_new_reference_available'] = False
            values = summary['components']['composite_R_model_kcal_mol']; value = values['median']
            missing = values['missing_members']; spread = values['range']
        else:
            if parent['pool'] != 'balanced': raise InvalidArtifact('unknown archived pool')
            arms = [next(p for p in pools if p['root_case_id'] == parent['root_case_id'] and p['pool'] == name) for name in ('La4', 'Ca5')]
            vals = [p['methods']['context_union']['R'] for p in arms]
            value = None if any(v is None for v in vals) else sum(vals)/2
            missing = sorted(set(cid for p in arms for cid in p['methods']['context_union']['missing_members'])); spread = None
            summary = {'equal_mean_of_complete_La4_and_Ca5_medians': True}
        call = decision(value, ref['bands']); methods = {k: dict(v) for k, v in parent['methods'].items()}
        methods['context_union'] = {'R': value, 'decision': call, 'outcome': outcome(call, parent['expected_class']),
                                   'missing_members': missing, 'required_members': len(members), 'within_protein_range': spread}
        pools.append({**parent, 'methods': methods, 'union_component_summary': summary})
    for triple in old['triples']:
        if triple['root_case_id'] not in roots: continue
        members = [next(r for r in selected if r['case_id'] == cid) for cid in triple['members']]
        if len(members) != 3 or any(r['source_conditioning_metal'] != 'La' for r in members):
            raise InvalidArtifact('archived triple scope changed')
        summary = pool_summary(members, component_idx, pool_bands); v = summary['components']['composite_R_model_kcal_mol']
        call = decision(v['median'], ref['bands']); methods = {k: dict(value) for k, value in triple['methods'].items()}
        methods['context_union'] = {'R': v['median'], 'decision': call, 'outcome': outcome(call, triple['expected_class']),
                                   'missing_members': v['missing_members'], 'required_members': 3, 'within_protein_range': v['range']}
        triples.append({**triple, 'methods': methods})
    if len(pools) != 3*len(roots) or len(triples) != 4*len(roots): raise InvalidArtifact('strict-summary denominator changed')
    def count(rr, method):
        counts = Counter(r['methods'][method]['outcome'] for r in rr)
        return {'denominator': len(rr), **{k: counts[k] for k in ('correct', 'wrong', 'inconclusive', 'unavailable')}}
    methods = list(rows[0]['methods']); subsets = {'all': rows,
          'La_noncanonical': [r for r in rows if r['source_conditioning_metal'] == 'La'],
          'Ca': [r for r in rows if r['source_conditioning_metal'] == 'Ca']}
    counts = {'single_sources': {name: {m: count(rr, m) for m in methods} for name, rr in subsets.items()},
              'strict_pools': {name: {m: count([r for r in pools if r['pool'] == name], m) for m in methods} for name in ('La4', 'Ca5', 'balanced')},
              'La_triples': {m: count(triples, m) for m in methods}}
    matched = {}
    for baseline_method in ('context_composite', 'DFT', 'core_native'):
        common = [r for r in rows if r['methods'][baseline_method]['R'] is not None and r['methods']['context_union']['R'] is not None]
        matched[baseline_method] = {'baseline': count(common, baseline_method), 'union': count(common, 'context_union'),
                 'transitions': [r for r in common if r['methods'][baseline_method]['outcome'] != r['methods']['context_union']['outcome']]}
    final = {'protocol_id': PROTOCOL, 'scope': scope, 'collection': record(collection), 'reference': record(reference_path),
             'baseline': record(baseline), 'sources': record(sources), 'rows': rows, 'pools': pools, 'triples': triples,
             'counts': counts, 'matched': matched, 'group_denominator': len(roots), 'new_molecular_calls': 0,
             'interpretation': 'Consumed structural replicas; union changes preparation/composition/cavity. Strict medians are not thermal populations.',
             'analysis_implementation': record(__file__)}
    write_new(output, final)
    return {k: final[k] for k in ('counts', 'group_denominator')}


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='op', required=True)
    q = sub.add_parser('collect'); q.add_argument('--stage', required=True); q.add_argument('--output', required=True)
    q = sub.add_parser('reference'); q.add_argument('--collection', required=True); q.add_argument('--output', required=True)
    q = sub.add_parser('report')
    for k in ('collection', 'reference-path', 'static-reference', 'output'): q.add_argument('--'+k, required=True)
    q = sub.add_parser('transfer-compare')
    for k in ('collection', 'reference-path', 'baseline', 'sources', 'output'): q.add_argument('--'+k, required=True)
    q.add_argument('--scope', choices=('pilot36', 'primary225'), default='pilot36')
    args = vars(p.parse_args()); op = args.pop('op').replace('-', '_'); print(json.dumps(globals()[op](**args), indent=2))


if __name__ == '__main__': main()
