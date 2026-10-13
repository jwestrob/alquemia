#!/usr/bin/env python3
"""Merge the two actual alternative-source batches without changing old records."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import time
import traceback

HERE = Path(__file__).resolve().parent
POLICIES = ('operational', 'mathematical', 'origin')


def pin(path):
    path = Path(path).resolve(); h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''): h.update(block)
    return {'path': str(path), 'sha256': h.hexdigest()}


def read(value):
    if isinstance(value, dict):
        assert pin(value['path']) == value, value['path']; value = value['path']
    return json.loads(Path(value).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def freeze(args):
    batches = [read(path) for path in args.batch_config]
    assert sorted(b['denominator'] for b in batches) == [5, 10]
    for key in ('scope', 'primary', 'prior_overlay', 'canonical_reference'):
        assert batches[0][key] == batches[1][key]
    ids = [sid for b in batches for sid in b['declared_sources']]
    assert len(set(ids)) == len(ids) == 15
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=False)
    config = {'protocol_id': 'Nikasha_two_Ca_full15_alternative_source_merge_v1',
        'implementation': pin(__file__), 'wrapper': pin(HERE/'merge.sbatch'),
        **{key: batches[0][key] for key in ('scope', 'primary', 'prior_overlay', 'canonical_reference')},
        'batches': [{'config': pin(path), 'expected_overlay': str(Path(batch['output'])/'FINAL_COLLECTION_ALTERNATIVE_SOURCE_OVERLAY.json')}
                    for path, batch in zip(args.batch_config, batches)],
        'declared_new_sources': ids, 'output': str(out), 'new_fit': False,
        'new_triples': False, 'new_molecular_calls': 0, 'production_changed': False}
    write(out/'manifest.json', config); print(json.dumps({'manifest': pin(out/'manifest.json'), 'executed': False}))


def run(path):
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('Actual Slurm required for scientific collection/report')
    config = read(path); assert config['implementation'] == pin(__file__)
    assert config['wrapper'] == pin(HERE/'merge.sbatch')
    out = Path(config['output']); started = time.monotonic()
    receipt = {'manifest': pin(path), 'job_id': os.environ['SLURM_JOB_ID'], 'status': 'started', 'new_molecular_calls': 0}
    try:
        primary = read(config['primary']); prior = read(config['prior_overlay']); scope = read(config['scope'])
        assert len(primary['rows']) == len(prior['rows']) == len(scope['rows']) == 125
        assert prior['primary'] == config['primary'] and [r['primary'] for r in prior['rows']] == primary['rows']
        reference = primary['reference']
        assert reference == prior['reference'] == scope['reference']
        assert config['canonical_reference'] == reference['source']
        assert read(config['canonical_reference'])['calibration'] == reference['frozen_calibration']
        replacements = {}; batch_records = []; scientific_costs = []
        for batch in config['batches']:
            b = read(batch['config']); overlay = read(batch['expected_overlay'])
            assert b['primary'] == config['primary'] and b['prior_overlay'] == config['prior_overlay']
            assert overlay['primary'] == config['primary'] and overlay['prior_recovery_overlay'] == config['prior_overlay']
            assert overlay['reference'] == reference and len(overlay['rows']) == 125
            assert [r['primary'] for r in overlay['rows']] == primary['rows']
            actual = read(overlay['two_Ca_collection'])
            assert actual['reference'] == reference and not actual['new_calibration'] and not actual['new_Ca_triples']
            by_new = {r['source_case_id']: r for r in actual['rows']}
            assert set(by_new) == set(b['declared_sources']) and len(by_new) == b['denominator']
            found = []
            for old, current in zip(prior['rows'], overlay['rows']):
                assert {k: v for k, v in current.items() if k != 'two_Ca_source_recovery'} == old
                replacement = current['two_Ca_source_recovery']
                if replacement is None: continue
                original_id = replacement['original_source_id']; new_id = replacement['source_case_id']
                assert original_id == old['primary']['source_case_id'] and original_id not in replacements
                assert old['primary']['decisions']['operational']['outcome'] == 'unavailable'
                assert old.get('technical_recovery') is None and old.get('preparation_recovery') is None
                assert all(replacement[key] == value for key, value in by_new[new_id].items())
                assert replacement['expected_class'] == old['primary']['expected_class']
                replacements[original_id] = replacement; found.append(new_id)
            assert set(found) == set(b['declared_sources'])
            native = read(b['native_collection'])
            assert native['manifest'] == b['native_master'] and native['denominator'] == b['denominator']
            batch_records.append({'config': batch['config'], 'overlay': pin(batch['expected_overlay']),
                                  'collection': overlay['two_Ca_collection'], 'native_collection': pin(b['native_collection'])})
            scientific_costs.append({'denominator': b['denominator'], 'native_lanes': native['lanes'],
                                     'solvent_backend_cost': actual['cost'],
                                     'cost_scope': 'Actual two-Ca scoring batch only; excludes folding, preparation and earlier campaigns.'})
        assert len(replacements) == 15
        assert {r['source_case_id'] for r in replacements.values()} == set(config['declared_new_sources'])
        groups = Counter(r['biological_group'] for r in replacements.values())
        assert groups == {'q9z4j7': 5, 'c5axv8': 5, 'q88jh5': 5}
        rows = [{**old, 'two_Ca_source_recovery': replacements.get(old['primary']['source_case_id'])} for old in prior['rows']]
        selected = [r['two_Ca_source_recovery'] or r.get('preparation_recovery') or r.get('technical_recovery') or r['primary'] for r in rows]
        counts = {policy: dict(Counter(r['decisions'][policy]['outcome'] for r in selected)) for policy in POLICIES}
        common = [r for r in selected if all(r['decisions'][p]['outcome'] != 'unavailable' for p in ('operational', 'origin'))]
        result = {'protocol_id': config['protocol_id'], 'manifest': pin(path), 'primary': config['primary'],
            'prior_recovery_overlay': config['prior_overlay'], 'batches': batch_records, 'rows': rows,
            'declared_sources': 125, 'new_two_Ca_sources': 15, 'new_sources_by_group': dict(groups),
            'primary_counts': primary['counts'], 'prior_recovery_counts': prior['annotated_recovery_counts'],
            'annotated_alternative_source_counts': counts,
            'common_origin_operational_coverage': {'sources': len(common), 'counts': {p: dict(Counter(r['decisions'][p]['outcome'] for r in common)) for p in ('operational', 'origin')}},
            'reference': reference, 'actual_scoring_costs': scientific_costs,
            'new_fit': False, 'new_triples': False, 'production_changed': False,
            'interpretation': 'Original125 unchanged; three prior technical/preparation recoveries and all15 new two-Ca source hypotheses are separately retained. Correlated samples are not independent proteins; no affinity or global spectator response is inferred.'}
        write(out/'RESULT.json', result)
        fields = ['original_source_id', 'selected_source_id', 'source_kind', 'biological_group', 'expected_class',
                  'origin_outcome', 'operational_outcome', 'mathematical_outcome', 'origin_R_kcal_mol',
                  'operational_R_kcal_mol', 'mathematical_R_kcal_mol', 'preparation_reason', 'execution_reason']
        with (out/'SOURCES.tsv').open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
            for item, selected_row in zip(rows, selected):
                kind = next((k for k in ('two_Ca_source_recovery', 'preparation_recovery', 'technical_recovery') if item.get(k) is not None), 'primary')
                record = {'original_source_id': item['primary']['source_case_id'], 'selected_source_id': selected_row['source_case_id'],
                    'source_kind': kind, 'biological_group': selected_row['biological_group'], 'expected_class': selected_row['expected_class'],
                    'preparation_reason': selected_row.get('preparation_reason'), 'execution_reason': selected_row.get('execution_reason')}
                for policy in POLICIES:
                    record[policy+'_outcome'] = selected_row['decisions'][policy]['outcome']
                    point = selected_row['pool'].get(policy)
                    record[policy+'_R_kcal_mol'] = point['R_model_kcal_mol'] if point is not None else None
                writer.writerow(record)
        lines = ['# All fifteen two-Ca alternative sources', '',
            'Both actual batches are included. Original125 records and the earlier three recoveries remain unchanged.', '',
            '| Policy | Correct | Wrong | Inconclusive | Unavailable |', '|---|---:|---:|---:|---:|']
        for policy in POLICIES:
            c = counts[policy]; lines.append('| '+policy+' | '+' | '.join(str(c.get(k, 0)) for k in ('correct','wrong','inconclusive','unavailable'))+' |')
        lines += ['', 'These are separately recorded alternative-source outcomes under the frozen calibration, not fifteen new independent proteins or a global spectator model.', '',
                  'SOURCES.tsv/RESULT.json retain all decisions, original records and actual batch costs. No new fit or triple policy was introduced.', '',
                  'Costs cover the actual native/solvent scoring batches only. Folding, preparation and earlier attempts are separate; no CPU/GPU speedup is inferred.']
        with (out/'REPORT.md').open('x') as handle: handle.write('\n'.join(lines)+'\n')
        receipt.update(status='complete', result=pin(out/'RESULT.json'), table=pin(out/'SOURCES.tsv'), report=pin(out/'REPORT.md'))
    except BaseException as exc:
        receipt.update(status='failed', reason=str(exc), traceback=traceback.format_exc()); raise
    finally:
        used = resource.getrusage(resource.RUSAGE_SELF)
        receipt.update(wall_seconds=time.monotonic()-started, process_CPU_seconds=used.ru_utime+used.ru_stime,
                       peak_RSS_KiB=used.ru_maxrss)
        write(out/('EXECUTION_'+os.environ['SLURM_JOB_ID']+'.json'), receipt)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='op', required=True)
    p = sub.add_parser('freeze'); p.add_argument('--batch-config', action='append', required=True); p.add_argument('--output', required=True)
    p = sub.add_parser('run'); p.add_argument('--manifest', required=True)
    args = parser.parse_args()
    if args.op == 'freeze': freeze(args)
    else: print(json.dumps(run(args.manifest), indent=2))
