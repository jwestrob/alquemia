#!/usr/bin/env python3
"""Summarize an immutable completed comparison; no new scoring or calibration."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
import time


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def checked(p):
    assert pin(p['path']) == p, p
    return json.loads(Path(p['path']).read_text())


def table(path, rows):
    with path.open('x') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--result', required=True)
    p.add_argument('--sha', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Slurm allocation required')
    started = time.monotonic()
    source = pin(a.result)
    assert source['sha256'] == a.sha
    result = checked(source)
    assert result['status'] == 'complete' and result['stage'] == 'transfer'
    output = Path(a.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    lookup = {r['record_id']: r for r in result['rows']}
    cohorts = {'expanded135': [r for r in result['rows'] if r['family'] == 'expanded135'],
               'primary_ca125': [r for r in result['rows'] if r['family'] == 'primary_ca125'],
               'current_ca125_annotated_overlay': [lookup[r['record_id']] for r in result['current_ca125_annotated_overlay']],
               'two_ca15': [r for r in result['rows'] if r['family'] in ('two_ca_q9_5', 'two_ca_c5q88_10')]}
    flat, groups, changes = [], [], []
    for row in result['rows']:
        for policy in ('operational', 'mathematical', 'origin'):
            new, old = row['pool'][policy], row['baseline']['pool'][policy]
            record = {k: row[k] for k in ('record_id', 'case_id', 'source_case_id', 'protein', 'biological_group',
                                          'family', 'role', 'expected_class', 'source_conditioning_metal')}
            record.update(policy=policy,
                          candidate_R=None if new is None else new['R_model_kcal_mol'],
                          baseline_R=None if old is None else old['R_model_kcal_mol'],
                          candidate_decision=row['decisions'][policy]['decision'],
                          candidate_outcome=row['decisions'][policy]['outcome'],
                          baseline_decision=row['baseline']['decisions'][policy]['decision'],
                          baseline_outcome=row['baseline']['decisions'][policy]['outcome'],
                          old_band_candidate_outcome=row['old_band_transfer'][policy]['outcome'],
                          missing_cells=json.dumps(row['missing_cells'], sort_keys=True))
            for z in ('Ca', 'La'):
                for label, entry in (('candidate', new), ('baseline', old)):
                    record[label+'_'+z+'_geometry'] = None if entry is None else entry[z]['candidate']
                    record[label+'_'+z+'_work_kcal_mol'] = None if entry is None else entry[z]['work_kcal_mol']
            flat.append(record)
    table(output/'SOURCES.tsv', flat)
    for name, rows in cohorts.items():
        by_group = defaultdict(list)
        for row in rows:
            by_group[(row['biological_group'], row['source_conditioning_metal'])].append(row)
        for (protein, metal), members in by_group.items():
            for policy in ('operational', 'mathematical', 'origin'):
                common = [r for r in members if r['pool'][policy] is not None and r['baseline']['pool'][policy] is not None]
                record = {'cohort': name, 'biological_group': protein, 'source_conditioning_metal': metal,
                          'policy': policy, 'declared': len(members), 'common_coverage': len(common),
                          'candidate_outcomes': dict(Counter(r['decisions'][policy]['outcome'] for r in members)),
                          'baseline_outcomes': dict(Counter(r['baseline']['decisions'][policy]['outcome'] for r in members))}
                for label in ('candidate', 'baseline'):
                    values = [(r if label == 'candidate' else r['baseline'])['pool'][policy]['R_model_kcal_mol'] for r in common]
                    record[label+'_common_min'] = min(values) if values else None
                    record[label+'_common_max'] = max(values) if values else None
                    record[label+'_common_range'] = max(values)-min(values) if values else None
                groups.append(record)
        for r in rows:
            if r['decisions']['operational'] != r['baseline']['decisions']['operational']:
                changes.append({'cohort': name, 'record_id': r['record_id'], 'protein': r['protein'],
                                'old': r['baseline']['decisions']['operational'], 'new': r['decisions']['operational']})
    native = [checked(result[k]) for k in ('canonical_native_collection', 'native_collection')]
    details, child_receipts = [], []
    for collection in native:
        for t in collection['tasks'].values():
            detail = checked(t['result'])
            details.append(detail)
            if not detail.get('reuse'):
                for medium in ('water', 'vacuum'):
                    if medium in detail:
                        child_receipts.append(checked(detail[medium]['execution']))
    residuals = [abs(c['residual_Eh']) for d in details for c in d.get('checks', [])]
    summary = {'result': source, 'implementation': pin(__file__), 'cohorts': result['cohorts'],
               'groups': groups, 'operational_transitions': changes,
               'triple_outcomes': {policy: dict(Counter(t['decisions'][policy]['outcome'] for t in result['triples']))
                                   for policy in ('operational', 'mathematical', 'origin')},
               'native': {'declared': len(details), 'available': sum(d['status'] == 'available' for d in details),
                          'unavailable': sum(d['status'] != 'available' for d in details),
                          'reused_inputs': sum(bool(d.get('reuse')) for d in details),
                          'new_component_calls': sum(c['native_component_calls'] for c in native),
                          'max_water_component_abs_residual_Eh': max(residuals) if residuals else None,
                          'sum_measured_child_CPU_seconds': sum(c.get('child_CPU_seconds', 0) for c in child_receipts),
                          'max_per_child_peak_RSS_KiB': max((c.get('peak_RSS_KiB', 0) for c in child_receipts), default=None)},
               'limits': ['Group ranges use common coverage; this is not a new group classification rule.',
                          'Cohorts overlap and structures are correlated; no fresh blind test.',
                          'Component replay timings exclude previously computed charges and MACE energies.',
                          'Measured child CPU and per-child RSS are not whole-allocation utilization or peak memory.'],
               'job_id': os.environ['SLURM_JOB_ID'], 'wall_seconds': time.monotonic()-started,
               'production_changed': False}
    with (output/'SUMMARY.json').open('x') as f:
        json.dump(summary, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({'native': summary['native'], 'triple_outcomes': summary['triple_outcomes'],
                      'operational_transitions': changes}, indent=2))


if __name__ == '__main__':
    main()
