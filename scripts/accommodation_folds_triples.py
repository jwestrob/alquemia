"""All prespecified three-of-four La-source subsets; consumed development only."""
from __future__ import annotations
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import statistics

from affordable_common import InvalidArtifact, read_json, record, write_new
from accommodation_folds_compare import decision

SCHEMA = 'PQQ_all_three_of_four_La_source_development_v1'
METHOD_FIELDS = {'native': 'native_R_model_kcal_mol', 'composite': 'composite_R_model_kcal_mol'}


def outcome(call, expected):
    if call == expected + '-supported': return 'correct'
    if call in ('Ca-supported', 'La-supported'): return 'wrong'
    return call


def summarize(comparison, plan, output):
    source = read_json(comparison); rows = []; proteins = []
    if source['schema_version'] != 'PQQ_fold_conditioning_comparison_v1' or source['protein_denominator'] != 25:
        raise InvalidArtifact('requires the completed declared 25-protein fold comparison')
    index = {(r['case_id'], r['representation']): r for r in source['rows']}
    if len(index) != 500: raise InvalidArtifact('original 250-source/two-representation scope changed')
    for parent in source['protein_summaries']:
        rep = parent['representation']; members = parent['pools']['La4']['declared_members']
        if len(members) != 4 or len(set(members)) != 4:
            raise InvalidArtifact('exact four original noncanonical La sources required')
        for cid in members:
            row = index[cid, rep]
            if row['canonical_coordinate_match'] or row['source_conditioning_metal'] != 'La' or row['root_case_id'] != parent['root_case_id']:
                raise InvalidArtifact('triple pool is not the original declared La4 pool')
        triples = []
        for selected in combinations(members, 3):
            result = {'root_case_id': parent['root_case_id'], 'representation': rep,
                      'expected_class': parent['expected_class'], 'members': list(selected),
                      'omitted_member': next(c for c in members if c not in selected), 'methods': {}}
            for method, field in METHOD_FIELDS.items():
                vals = [index[cid, rep][field] for cid in selected]
                missing = [cid for cid, value in zip(selected, vals) if value is None]
                value = None if missing else statistics.median(vals)
                call = decision(value, source['frozen_bands'][rep][method])
                result['methods'][method] = {'median_R_model_kcal_mol': value, 'decision': call,
                    'outcome': outcome(call, parent['expected_class']), 'missing_members': missing,
                    'required_members': 3, 'available_members': 3 - len(missing)}
            triples.append(result); rows.append(result)
        protein = {'root_case_id': parent['root_case_id'], 'representation': rep,
                   'expected_class': parent['expected_class'], 'declared_La4_members': members,
                   'triples': triples, 'methods': {}}
        for method in METHOD_FIELDS:
            cc = Counter(t['methods'][method]['outcome'] for t in triples)
            counts = {k: cc[k] for k in ('correct', 'wrong', 'inconclusive', 'unavailable')}
            worst = next(k for k in ('wrong', 'unavailable', 'inconclusive', 'correct') if counts[k])
            scored_worst = next((k for k in ('wrong', 'inconclusive', 'correct') if counts[k]), 'unavailable')
            protein['methods'][method] = {'counts': counts, 'denominator': 4,
                'all_four_correct': counts['correct'] == 4, 'any_wrong': counts['wrong'] > 0,
                'any_inconclusive': counts['inconclusive'] > 0, 'any_unavailable': counts['unavailable'] > 0,
                'worst_outcome': worst, 'worst_scored_outcome': scored_worst}
        proteins.append(protein)
    if len(rows) != 200 or len(proteins) != 50:
        raise InvalidArtifact('three-fold descriptor denominator changed')
    totals = {}
    for rep in ('core', 'context'):
        totals[rep] = {}
        for method in METHOD_FIELDS:
            calls = Counter(r['methods'][method]['outcome'] for r in rows if r['representation'] == rep)
            pp = [p['methods'][method] for p in proteins if p['representation'] == rep]
            totals[rep][method] = {'triple_denominator': 100, 'protein_denominator': 25,
                'triples': {k: calls[k] for k in ('correct', 'wrong', 'inconclusive', 'unavailable')},
                'proteins_all_four_correct': sum(p['all_four_correct'] for p in pp),
                'proteins_any_wrong': sum(p['any_wrong'] for p in pp),
                'proteins_any_inconclusive': sum(p['any_inconclusive'] for p in pp),
                'proteins_any_unavailable': sum(p['any_unavailable'] for p in pp),
                'protein_worst_outcomes': dict(Counter(p['worst_outcome'] for p in pp))}
    result = {'schema_version': SCHEMA, 'comparison': record(comparison), 'plan': record(plan),
              'implementation': record(__file__), 'frozen_bands': source['frozen_bands'],
              'rows': rows, 'proteins': proteins, 'counts': totals,
              'new_molecular_calls': 0, 'threshold_refitted': False, 'production_changed': False,
              'new_PLM_validation': False, 'all_evidence_consumed': True,
              'interpretation': 'every three-of-four subset on consumed reference structures; not independent biological or PLM validation'}
    write_new(output, result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('comparison', 'plan', 'output'): p.add_argument('--' + name, required=True)
    result = summarize(**vars(p.parse_args())); print(json.dumps(result['counts']))


if __name__ == '__main__': main()
