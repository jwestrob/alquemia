"""Join the completed, frozen-reference comparisons; perform no molecular work."""
from __future__ import annotations
import argparse
from collections import Counter
from itertools import combinations
from pathlib import Path

from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from accommodation_fold_proposals import summarize_rows
from consistent_context import PROTOCOL as UNION_PROTOCOL
from nikasha_pool import SCALED_PROTOCOL

METHODS = ('context_composite', 'DFT', 'context_union', 'adaptive_pool')
LABELS = {'context_composite': 'Released MACE/GFN2', 'DFT': 'Preserved DFT',
          'context_union': 'Consistent pocket', 'adaptive_pool': 'Adaptive common pool'}


def key(row, kind):
    if kind == 'rows': return row['case_id']
    return (row['root_case_id'], row.get('pool'), tuple(sorted(row['members'])))


def join_rows(base, extra, kind, method, target):
    """Retain baseline values and reject different populations or biological labels."""
    index = {key(r, kind): r for r in extra}
    if len(index) != len(extra) or set(index) != {key(r, kind) for r in base}:
        raise InvalidArtifact('comparison populations differ: ' + kind)
    result = []
    for original in base:
        other = index[key(original, kind)]
        if (original['expected_class'] != other['expected_class'] or
                original['root_case_id'] != other['root_case_id']):
            raise InvalidArtifact('source group or label changed')
        if kind == 'rows' and original['source_conditioning_metal'] != other['source_conditioning_metal']:
            raise InvalidArtifact('source conditioning changed')
        for name in ('context_composite', 'DFT'):
            if any(original['methods'][name][f] != other['methods'][name][f]
                   for f in ('R', 'decision', 'outcome')):
                raise InvalidArtifact('archived baseline result changed')
        result.append({**original, 'methods': {**original['methods'], target: other['methods'][method]}})
    return result


def compare(context, adaptive, output):
    union, pool = read_json(context), read_json(adaptive)
    if (union['protocol_id'] != UNION_PROTOCOL or union['scope'] != 'primary225' or
            pool['protocol_id'] != SCALED_PROTOCOL or union['baseline'] != pool['prior_comparison'] or
            union['sources'] != pool['sources']):
        raise InvalidArtifact('different source panel, baseline or candidate protocols')
    base = read_json(verify(union['baseline']))
    ur, ar = read_json(verify(union['reference'])), read_json(verify(pool['reference']))
    if ur['status'] != 'available' or any(v['status'] != 'available' for v in ar['variants'].values()):
        raise InvalidArtifact('candidate decision reference unavailable')
    merged = {}
    for kind, expected in (('rows', 225), ('pools', 75), ('triples', 100)):
        if len(base[kind]) != expected: raise InvalidArtifact('original denominator changed')
        rows = join_rows(base[kind], union[kind], kind, 'context_union', 'context_union')
        merged[kind] = join_rows(rows, pool[kind], kind, 'pool_operational_new', 'adaptive_pool')
    rows, pools, triples = (merged[k] for k in ('rows', 'pools', 'triples'))
    subsets = {'all225': rows,
               'La100': [r for r in rows if r['source_conditioning_metal'] == 'La'],
               'Ca125': [r for r in rows if r['source_conditioning_metal'] == 'Ca'],
               **{name: [r for r in pools if r['pool'] == name] for name in ('La4', 'Ca5', 'balanced')},
               'La100_triples': triples}
    counts = {name: summarize_rows(group, METHODS) for name, group in subsets.items()}
    matched = {}
    for name, group in subsets.items():
        matched[name] = {}
        for left, right in combinations(METHODS, 2):
            common = [r for r in group if all(r['methods'][m]['R'] is not None and
                      r['methods'][m]['decision'] != 'unavailable' for m in (left, right))]
            matched[name][left + '__' + right] = {
                'declared': len(group), 'common': len(common),
                'counts': summarize_rows(common, (left, right)),
                'outcome_transitions': dict(Counter(r['methods'][left]['outcome'] + '->' +
                    r['methods'][right]['outcome'] for r in common))}
    transitions = []
    for r in rows:
        for method in ('context_union', 'adaptive_pool'):
            before, after = r['methods']['context_composite'], r['methods'][method]
            if before['outcome'] != after['outcome']:
                transitions.append({'case_id': r['case_id'], 'root_case_id': r['root_case_id'],
                    'source_conditioning_metal': r['source_conditioning_metal'],
                    'expected_class': r['expected_class'], 'candidate': method,
                    'released': before, 'new': after})
    result = {'context_comparison': record(context), 'adaptive_comparison': record(adaptive),
              'baseline': union['baseline'], 'sources': union['sources'],
              'references': {'context_union': union['reference'], 'adaptive_pool': pool['reference']},
              'method_labels': LABELS, **merged, 'counts': counts, 'matched': matched,
              'source_transitions': transitions, 'implementation': record(__file__),
              'new_molecular_calls': 0, 'new_thresholds_fitted': False, 'production_changed': False,
              'interpretation': 'Consumed structures; fixed per-method canonical references. '
                                'Strict groups require all members. Correlated triples are not independent proteins.'}
    write_new(output, result)
    lines = ['# Nikasha: completed next-phase comparison', '',
             'Each changed scorer uses its own frozen canonical-only reference. '
             'Counts retain all failures; abstention is separate from a correct decision.', '',
             '| Population | Method | Correct | Wrong | Inconclusive | Unavailable | Total |',
             '|---|---|---:|---:|---:|---:|---:|']
    for name, tally in counts.items():
        for method in METHODS:
            c = tally[method]
            lines.append('| ' + name + ' | ' + LABELS[method] + ' | ' +
                         ' | '.join(str(c[k]) for k in ('correct', 'wrong', 'inconclusive', 'unavailable', 'denominator')) + ' |')
    lines += ['', 'The JSON retains every source, strict group, correlated triple, changed outcome, '
              'and pairwise common-coverage comparison. Raw scales are protocol-specific; '
              'these outputs are not binding free energies or physiological occupancy predictions.', '']
    with Path(output).with_suffix('.md').open('x') as f: f.write('\n'.join(lines))
    return counts['all225']


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('context', 'adaptive', 'output'): p.add_argument('--' + name, required=True)
    print(compare(**vars(p.parse_args())))
