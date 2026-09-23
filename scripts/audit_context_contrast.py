"""Decompose saved three-source versus ten-fold contrasts; no molecular work."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
from statistics import median

from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new
from mace_hybrid import EV_TO_KCAL


def origin(matrix):
    ca = matrix['Ca']['origin']['components']
    la = matrix['La']['origin']['components']
    native = (ca['MACE_eV'] - la['MACE_eV']) * EV_TO_KCAL
    solvent = ((ca['GFN2_ALPB_hartree'] - ca['GFN2_vacuum_hartree']) -
               (la['GFN2_ALPB_hartree'] - la['GFN2_vacuum_hartree'])) * HA_TO_KCAL
    return {'native': native, 'solvent': solvent, 'total': native + solvent}


def selected(pool):
    p = pool['operational']
    result = dict(native=p['native_R_model_kcal_mol'],
                  solvent=p['solvation_delta_R_kcal_mol'],
                  total=p['composite_R_model_kcal_mol'])
    if abs(result['native'] + result['solvent'] - result['total']) > 1e-6:
        raise InvalidArtifact('saved selected-component closure differs')
    return result


def audit(triples, tenfold, output):
    new = read_json(triples); old = read_json(tenfold)
    if new['pair_denominator'] != 125 or len(new['rows']) != 125:
        raise InvalidArtifact('declared complete three-source transfer required')
    old_rows = {r['case_id']: r for r in old['rows']}
    if len(old_rows) != 225:
        raise InvalidArtifact('declared complete ten-fold transfer required')
    changed = [r for r in new['rows'] if not r['pool_reused']]
    if len(changed) != 55:
        raise InvalidArtifact('expected all55 changed source/context pairs')
    collections = {}; rows = []
    for r in changed:
        previous = old_rows[r['case_id']]
        pin = previous['precision_collection']; path = verify(pin)
        if str(path) not in collections:
            collections[str(path)] = {c['case_id']: c for c in read_json(path)['cases']}
        prior = collections[str(path)][r['case_id']]
        before = origin(prior['matrix']); after = origin(r['matrix'])
        before_pool = selected(prior['pool']); after_pool = selected(r['pool'])
        if abs(before_pool['total'] - previous['methods']['union_precision']['R']) > 1e-6:
            raise InvalidArtifact('ten-fold collection and comparison disagree')
        if abs(after_pool['total'] - r['methods']['triple_union']['R']) > 1e-6:
            raise InvalidArtifact('three-source collection and comparison disagree')
        row = {k: r[k] for k in ('pair_id', 'case_id', 'root_case_id', 'expected_class',
                                  'atom_count', 'tenfold_atom_count')}
        row.update(old_outcome=previous['methods']['union_precision']['outcome'],
                   new_outcome=r['methods']['triple_union']['outcome'])
        for component in ('native', 'solvent', 'total'):
            row['origin_' + component + '_delta'] = after[component] - before[component]
            row['pool_' + component + '_delta'] = after_pool[component] - before_pool[component]
            row['accommodation_' + component + '_delta'] = (
                row['pool_' + component + '_delta'] - row['origin_' + component + '_delta'])
        rows.append(row)
    metrics = {}
    for prefix in ('origin', 'pool', 'accommodation'):
        metrics[prefix] = {
            'median_absolute_delta_kcal_mol': {
                c: median(abs(r[f'{prefix}_{c}_delta']) for r in rows)
                for c in ('native', 'solvent', 'total')},
            'solvent_absolute_change_exceeds_native_count': sum(
                abs(r[f'{prefix}_solvent_delta']) > abs(r[f'{prefix}_native_delta']) for r in rows)}
    result = {'protocol_id': 'Nikasha_saved_context_contrast_decomposition_v1',
              'three_source_comparison': record(triples), 'tenfold_comparison': record(tenfold),
              'source_context_pairs': len(rows), 'distinct_structures': len({r['case_id'] for r in rows}),
              'protein_groups': len({r['root_case_id'] for r in rows}), 'units': 'model kcal/mol',
              'direction': 'three-source minus ten-fold; operational selection; no new reference',
              'interpretation': 'Descriptive correlated development pairs; no unique physical cause or statistical test.',
              'new_molecular_calls': 0, 'metrics': metrics, 'rows': rows,
              'implementation': record(__file__)}
    out = Path(output); out.mkdir(parents=True, exist_ok=False)
    write_new(out/'RESULT.json', result)
    with (out/'rows.csv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    return {k: v for k, v in result.items() if k != 'rows'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--triples', required=True)
    parser.add_argument('--tenfold', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.triples, args.tenfold, args.output), indent=2))
