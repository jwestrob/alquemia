"""Describe frozen-band transfer using archived scores; no new model calls."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import resource
import statistics
import time
import traceback


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def write(path, data):
    with Path(path).open('x') as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write('\n')


def describe(values):
    return {'n': len(values), 'minimum': min(values), 'median': statistics.median(values),
            'maximum': max(values)} if values else {'n': 0}


def run(args):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Scientific scalar analysis requires Slurm')
    start = time.monotonic()
    config = json.loads(args.config.read_text())
    out = Path(config['output'])
    out.mkdir(parents=True, exist_ok=False)
    receipt = {'job_id': os.environ['SLURM_JOB_ID'], 'new_molecular_calls': 0,
               'config': pin(args.config), 'status': 'running'}
    try:
        for item in config['pins']:
            assert pin(item['path']) == item, item['path']
        final = json.loads(Path(config['final']).read_text())
        canonical = json.loads(Path(config['canonical']).read_text())
        bands = final['reference']['frozen_calibration']
        assert final['reference']['source'] == pin(config['canonical'])
        rows = []
        for wrapper in final['rows']:
            r = (wrapper.get('two_Ca_source_recovery') or wrapper.get('preparation_recovery')
                 or wrapper.get('technical_recovery') or wrapper['primary'])
            rows.append(('recovered125', r))
        canonical_rows = [r for r in canonical['rows'] if r.get('role') == 'canonical_calibration']
        assert len(canonical_rows) == 25
        rows.extend(('canonical25', r) for r in canonical_rows)
        normalized = []
        for scope, r in rows:
            item = {'scope': scope, 'source_id': r['source_case_id'],
                    'group': r['biological_group'], 'expected_class': r['expected_class'],
                    'conditioning': r.get('source_conditioning_metal'), 'policies': {}}
            for policy in ('origin', 'operational'):
                decision = r['decisions'][policy]
                value = ((r.get('pool') or {}).get(policy) or {}).get('R_model_kcal_mol')
                b = bands[policy]['bands']
                row = {'R_kcal_mol': value, 'outcome': decision['outcome']}
                if value is not None:
                    ca_margin = b['Ca_max'] - value
                    la_margin = value - b['La_min']
                    row.update(Ca_band_margin_kcal_mol=ca_margin,
                               La_band_margin_kcal_mol=la_margin,
                               expected_band_margin_kcal_mol=ca_margin if r['expected_class'] == 'Ca' else la_margin)
                    inferred = 'Ca' if ca_margin >= 0 else 'La' if la_margin >= 0 else 'inconclusive'
                    outcome = 'inconclusive' if inferred == 'inconclusive' else 'correct' if inferred == r['expected_class'] else 'wrong'
                    assert outcome == decision['outcome'], (item['source_id'], policy)
                item['policies'][policy] = row
            a, b = item['policies']['origin'], item['policies']['operational']
            if a['R_kcal_mol'] is not None and b['R_kcal_mol'] is not None:
                item['accommodation_delta_R_kcal_mol'] = b['R_kcal_mol'] - a['R_kcal_mol']
                item['expected_band_margin_change_kcal_mol'] = b['expected_band_margin_kcal_mol'] - a['expected_band_margin_kcal_mol']
            normalized.append(item)
        summaries = {}
        for scope in ('canonical25', 'recovered125'):
            rr = [r for r in normalized if r['scope'] == scope]
            summaries[scope] = {'sources': len(rr), 'groups': len({r['group'] for r in rr}),
                'transitions': dict(Counter(r['policies']['origin']['outcome'] + '->' + r['policies']['operational']['outcome'] for r in rr)),
                'policy_counts': {p: dict(Counter(r['policies'][p]['outcome'] for r in rr)) for p in ('origin', 'operational')},
                'by_expected_class': {}}
            for cls in ('Ca', 'La'):
                members = [r for r in rr if r['expected_class'] == cls]
                groups = defaultdict(list)
                for r in members:
                    if 'accommodation_delta_R_kcal_mol' in r:
                        groups[r['group']].append(r['accommodation_delta_R_kcal_mol'])
                summaries[scope]['by_expected_class'][cls] = {
                    'sources': len(members), 'groups': len(groups),
                    'delta_R_by_source': describe([v for vv in groups.values() for v in vv]),
                    'delta_R_group_medians': {g: statistics.median(v) for g, v in groups.items()},
                    'expected_band_margin_change_by_source': describe([r['expected_band_margin_change_kcal_mol'] for r in members if 'expected_band_margin_change_kcal_mol' in r])}
        result = {'protocol': 'Nikasha_PedE_frozen_band_diagnosis_v1', 'pins': config['pins'],
                  'bands': bands, 'band_edge_shifts': {key: bands['operational']['bands'][key] - bands['origin']['bands'][key] for key in ('Ca_max', 'La_min')},
                  'summaries': summaries, 'rows': normalized, 'new_fit': False,
                  'interpretation': 'Descriptive consumed-source analysis; margins are distances to frozen decision bands, not probabilities or affinities.'}
        write(out/'RESULT.json', result)
        lines = ['# Frozen-band diagnosis', '', '| Scope | Source | Origin call | Accommodated call | Delta R | Change in expected-class band margin |',
                 '|---|---|---|---|---:|---:|']
        for r in normalized:
            if r['group'] == 'q88jh5':
                lines.append(f"| {r['scope']} | {r['source_id']} | {r['policies']['origin']['outcome']} | {r['policies']['operational']['outcome']} | {r.get('accommodation_delta_R_kcal_mol', float('nan')):.6f} | {r.get('expected_band_margin_change_kcal_mol', float('nan')):.6f} |")
        lines += ['', 'Units kcal/mol; positive expected-class margin is on the correct side of its frozen band.',
                  'This is descriptive and introduces no new bands, decisions, probabilities or independent biological observations.', '']
        with (out/'REPORT.md').open('x') as handle:
            handle.write('\n'.join(lines))
        receipt.update(status='complete', result=pin(out/'RESULT.json'), report=pin(out/'REPORT.md'))
    except Exception as error:
        receipt.update(status='failed', error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        receipt.update(elapsed_seconds=time.monotonic() - start,
                       process_CPU_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_stime,
                       peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(out/'RECEIPT.json', receipt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, required=True)
    run(parser.parse_args())
