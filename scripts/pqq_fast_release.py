"""Finite source-to-score fast PQQ release checks; no calibration fitting."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from affordable_common import InvalidArtifact, read_json, record, verify, write_new


def prepare_score(preparation, comparison, agreement, output, cpu_python, gpu_python):
    from compact_solvation_scanner import prepare_pairs
    p = read_json(preparation)
    if p['supported'] != p['denominator'] or any(c['status'] != 'prepared' for c in p['cases']):
        raise InvalidArtifact('source preparation has failures/mismatches; scores not silently reused')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    cp = read_json(verify(p['config']['calibration_implementation_pins']))
    pairs = {'schema_version': 'prepared_complete_context_pairs_v1', 'model': p['config']['model'],
             'software': p['config']['software'], 'orca': cp['orca_runtime']['executable'],
             'cases': p['cases'], 'source_preparation': record(preparation)}
    path = out / 'pairs.json'; write_new(path, pairs)
    return prepare_pairs(path, comparison, agreement, out / 'scoring', cpu_python, gpu_python)


def collect(preparation, result, output):
    p = read_json(preparation); r = read_json(result)
    old = read_json(verify(read_json(verify(r['manifest']))['archived_comparison']))
    index = {x['case_id']: x for x in r['rows']}; rows = []
    for case in p['cases']:
        cid = case['case_id']; row = index.get(cid)
        expected = case['expected_class'] + '-supported'
        previous = next(x for x in old['rows'] if x['case_id'] == cid and x['representation'] == 'context')
        rows.append({'case_id': cid, 'role': case['role'], 'biological_group': case['biological_group'],
                     'expected_class': case['expected_class'], 'preparation_status': case['status'],
                     'preparation_seconds': case['source_preparation_seconds'],
                     'source_reused': case['preparation_reused'], 'context_comparison': case['context_comparison'],
                     'score_status': row['status'] if row else 'missing',
                     'decision': row['published_PQQ_decision'] if row else None,
                     'classification_pass': bool(row and row['status'] == 'available' and row['published_PQQ_decision'] == expected),
                     'R_model_kcal_mol': row['composite_R_model_kcal_mol'] if row else None,
                     'archived_R_model_kcal_mol': previous['composite_R_model_kcal_mol'],
                     'score_difference_model_kcal_mol': row['composite_R_model_kcal_mol'] - previous['composite_R_model_kcal_mol'] if row and row['status'] == 'available' else None,
                     'endpoint_path_seconds': row['timing']['fresh_endpoint_path_seconds'] if row and row.get('timing') else None})
    result = {'protocol_id': r['protocol_id'], 'source_preparation': record(preparation), 'score_collection': record(result),
              'rows': rows, 'denominator': 28, 'calibration_denominator': 25, 'transfer_denominator': 3,
              'calibration_correct': sum(x['classification_pass'] for x in rows if x['role'] == 'calibration'),
              'transfer_correct': sum(x['classification_pass'] for x in rows if x['role'] != 'calibration'),
              'preparation_supported': p['supported'], 'source_preparation_seconds': p['fresh_preparation_seconds'],
              'source_reused': False, 'energy_reused': False, 'frozen_bands': r['published_context_bands'],
              'all_evidence_consumed': True, 'threshold_refitted': False, 'baseline_changed': False,
              'promotion_gate': len(rows) == 28 and p['supported'] == 28 and all(x['classification_pass'] for x in rows)}
    write_new(output, result)
    return {k: v for k, v in result.items() if k not in ('rows', 'frozen_bands')}


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='op', required=True)
    q = sub.add_parser('prepare-score')
    for name in ('preparation', 'comparison', 'agreement', 'output', 'cpu-python', 'gpu-python'):
        q.add_argument('--' + name, required=True)
    q = sub.add_parser('collect')
    for name in ('preparation', 'result', 'output'): q.add_argument('--' + name, required=True)
    a = vars(p.parse_args()); op = a.pop('op').replace('-', '_'); print(json.dumps(globals()[op](**a)))


if __name__ == '__main__': main()
