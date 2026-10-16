#!/usr/bin/env python3
"""Calibration-only rules applied unchanged to pinned, consumed PQQ scores."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time
import traceback

POLICIES = ('operational', 'mathematical')
RULES = ('existing_bands', 'canonical_midpoint', 'calibration_sensitivity_band')
OUTCOMES = ('correct', 'wrong', 'inconclusive', 'unavailable', 'calibration_unavailable')
IDENTITY = ('context_id', 'source_case_id', 'case_id', 'protein', 'biological_group',
            'source_conditioning_metal', 'expected_class')


def require(value, message):
    if not value:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def check(value):
    require(pin(value['path']) == value, 'Changed input: ' + value['path'])
    return Path(value['path'])


def read(value):
    return json.loads(check(value).read_text())


def write(path, data):
    with Path(path).open('x') as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write('\n')


def table(path, rows):
    with Path(path).open('x', newline='') as handle:
        fields = list(dict.fromkeys(k for row in rows for k in row))
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in row.items()})


def score(row, policy):
    cell = row['pools']['common6'][policy]
    value = None if cell is None else cell['R_model_kcal_mol']
    require(value is None or math.isfinite(value), 'Nonfinite saved score')
    return value


def fit(rows, policy):
    """This function receives calibration records only; no transfer argument."""
    require(len({r['record_id'] for r in rows}) == len(rows), 'Duplicate calibration record')
    ca = [score(r, policy) for r in rows if r['expected_class'] == 'Ca']
    la = [score(r, policy) for r in rows if r['expected_class'] == 'La']
    if not ca or not la or any(v is None for v in ca + la):
        return {'status': 'unavailable_calibration', 'midpoint': None}
    require(len(ca) + len(la) == len(rows), 'Unsupported calibration label')
    lower, upper = max(ca), min(la)
    if lower >= upper:
        return {'status': 'nonseparating_calibration', 'midpoint': None}
    return {'status': 'available', 'Ca_max': lower, 'La_min': upper,
            'midpoint': lower + (upper - lower) / 2,
            'training_record_ids': [r['record_id'] for r in rows],
            'Ca_boundary_record_ids': [r['record_id'] for r in rows if r['expected_class'] == 'Ca' and score(r, policy) == lower],
            'La_boundary_record_ids': [r['record_id'] for r in rows if r['expected_class'] == 'La' and score(r, policy) == upper]}


def classify(value, lower, upper, expected):
    if value is None:
        return {'decision': None, 'outcome': 'unavailable'}
    if lower is None or upper is None:
        return {'decision': None, 'outcome': 'calibration_unavailable'}
    call = 'Ca' if value < lower else 'La' if value > upper else 'inconclusive'
    return {'decision': call, 'outcome': 'inconclusive' if call == 'inconclusive' else 'correct' if call == expected else 'wrong'}


def counts(rows):
    values = Counter(r['outcome'] for r in rows)
    require(set(values) <= set(OUTCOMES), 'Unknown outcome')
    return {key: values[key] for key in OUTCOMES}


def compare(config, out):
    spec = importlib.util.spec_from_file_location('_old_decisions', check(config['old_scorer']))
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    inventory, canonical = read(config['inventory']), read(config['canonical'])
    designated = inventory['calibration_rule']['designated_record_ids']
    training = canonical['rows']
    require(len(designated) == 25 and [r['record_id'] for r in training] == designated, 'Canonical membership/order differs')
    inv = {r['record_id']: r for r in inventory['records']}
    for row in training:
        require(all(row[k] == inv[row['record_id']][k] for k in ('expected_class', 'biological_group', 'source_case_id', 'source_conditioning_metal')), 'Canonical identity/label differs')
    models, tests, held_out = {}, [], []
    for policy in POLICIES:
        full = fit(training, policy)
        frozen_bands = canonical['calibrations']['common6'][policy]['bands']
        require(full['status'] == 'available' and {k: full[k] for k in ('Ca_max', 'La_min')} == frozen_bands, 'Canonical extrema do not reproduce frozen bands')
        folds = []
        for group in sorted({r['biological_group'] for r in training}):
            retained = [r for r in training if r['biological_group'] != group]
            omitted = [r for r in training if r['biological_group'] == group]
            model = fit(retained, policy)
            require(not ({r['record_id'] for r in retained} & {r['record_id'] for r in omitted}), 'Held-out leakage')
            folds.append({'omitted_group': group, 'omitted_record_ids': [r['record_id'] for r in omitted], 'fit': model})
            for row in omitted:
                held_out.append({'policy': policy, 'group': group, 'record_id': row['record_id'],
                                 'score': score(row, policy), 'expected_class': row['expected_class'],
                                 'threshold': model['midpoint'], **classify(score(row, policy), model['midpoint'], model['midpoint'], row['expected_class'])})
        available = all(f['fit']['status'] == 'available' for f in folds)
        thresholds = [full['midpoint']] + [f['fit']['midpoint'] for f in folds] if available else []
        models[policy] = {'full': full, 'existing_bands': frozen_bands, 'leave_group_out': folds,
                          'sensitivity_bounds': [min(thresholds), max(thresholds)] if available else [None, None]}
        tests.append({'test': 'canonical_midpoint_inside_gap_and_exact_tie_abstains', 'policy': policy,
                      'pass': full['Ca_max'] < full['midpoint'] < full['La_min'] and classify(full['midpoint'], full['midpoint'], full['midpoint'], training[0]['expected_class'])['outcome'] == 'inconclusive'})
        incomplete = [r for r in training if r['expected_class'] == 'Ca']
        tests.append({'test': 'incomplete_real_calibration_copy_unavailable', 'policy': policy,
                      'pass': fit(incomplete, policy)['status'] == 'unavailable_calibration',
                      'fixture': 'Explicit subset of real calibration with all La rows removed; not new evidence'})
    # Write all thresholds before opening any transfer outcomes.
    write(out / 'FROZEN_THRESHOLDS.json', {'canonical': config['canonical'], 'designated_record_ids': designated, 'models': models})
    transfer = read(config['transfer'])
    require(transfer['declared_contexts'] == 257 and transfer['declared_memberships'] == 260, 'Transfer denominator changed')
    records = {r['record_id']: r for r in transfer['rows']}
    unique = {r['context_id']: r for r in transfer['unique_context_rows']}
    require(len(records) == 260 and len(unique) == 257, 'Duplicate or missing transfer records')
    canonical_contexts = {records[rid]['context_id'] for rid in designated}
    require(len(canonical_contexts) == 25, 'Canonical context mapping differs')
    for row in training:
        current = records[row['record_id']]
        require(all(row[k] == current[k] for k in ('expected_class', 'biological_group', 'case_id')), 'Canonical transfer identity differs')
        require(all(score(row, p) == score(current, p) for p in POLICIES), 'Canonical transfer score differs')
    for row in records.values():
        reference = unique[row['context_id']]
        require(row['record_id'] in reference['memberships'] and all(row[k] == reference[k] for k in IDENTITY), 'Alias identity differs')
        require(all(row['pools']['common6'][p] == reference['pools']['common6'][p] for p in POLICIES), 'Alias scientific result differs')

    def calls(value, policy, expected):
        model = models[policy]
        t = model['full']['midpoint']
        lower, upper = model['sensitivity_bounds']
        return {'existing_bands': old.decision(value, model['existing_bands'], expected),
                'canonical_midpoint': classify(value, t, t, expected),
                'calibration_sensitivity_band': classify(value, lower, upper, expected)}

    rows, strata = [], defaultdict(list)
    for cid, row in sorted(unique.items()):
        calibration = cid in canonical_contexts
        memberships = [records[rid] for rid in row['memberships']]
        for policy in POLICIES:
            value = score(row, policy)
            predictions = calls(value, policy, row['expected_class'])
            require(all(predictions['existing_bands'] == r['decisions_by_pool']['common6'][policy] for r in memberships), 'Old-band replay differs')
            for rule, decision in predictions.items():
                item = {k: row[k] for k in IDENTITY}
                item.update(policy=policy, rule=rule, score=value, calibration=calibration, record_ids=row['memberships'], **decision,
                            old_decision=predictions['existing_bands']['decision'], old_outcome=predictions['existing_bands']['outcome'])
                rows.append(item)
                groups = ['all_unique', 'canonical' if calibration else 'noncanonical',
                          'conditioning:' + str(row['source_conditioning_metal']), 'protein:' + row['biological_group']]
                for group in groups:
                    strata[group, policy, rule].append(item)
    by_context = {(r['context_id'], r['policy'], r['rule']): r for r in rows}
    ledger = read(transfer['historical_ledger'])
    cohorts = {'expanded135': ledger['cohorts']['expanded135'],
               'current_ca125': [r['record_id'] for r in ledger['cohorts']['current_ca125_annotated_overlay']],
               'two_ca15': ledger['cohorts']['two_ca15']}
    require({k: len(v) for k, v in cohorts.items()} == {'expanded135': 135, 'current_ca125': 125, 'two_ca15': 15}, 'Cohort denominators differ')
    require(set(records) == set(cohorts['expanded135']) | set(cohorts['current_ca125']), 'Cohort membership differs')
    require(len(ledger['records']) == 278 and len(transfer['superseded_primary_records']) == 18, 'Historical ledger differs')
    superseded = {r['record_id'] for r in transfer['superseded_primary_records']}
    require(len(superseded) == 18 and not (superseded & set(records)) and
            superseded == {r['record_id'] for r in ledger['records']} - set(records), 'Historical exclusions overlap current records')
    for cohort, ids in cohorts.items():
        for policy in POLICIES:
            for rule in RULES:
                strata['cohort_memberships:' + cohort, policy, rule] = [by_context[records[rid]['context_id'], policy, rule] for rid in ids]

    triples = []
    lookup = {r['case_id']: r for r in records.values() if r['family'] == 'expanded135'}
    declared_triples = inventory['cohorts']['existing_100_triples']
    require(len(declared_triples) == 100 and len(transfer['existing_triples']['common6']) == 100, 'Triple denominator changed')
    for declared, saved in zip(declared_triples, transfer['existing_triples']['common6']):
        require(all(declared[k] == saved[k] for k in ('triple_id', 'members', 'pool_case_ids')) and declared['status'] == saved['prior_preparation_status'], 'Triple definition changed')
        members = [lookup.get(k) for k in saved['pool_case_ids']]
        expected = next((r['expected_class'] for r in members if r is not None), None)
        require(all(r['expected_class'] == expected for r in members if r is not None), 'Mixed triple labels')
        for policy in POLICIES:
            available = declared['status'] == 'prepared' and all(r is not None and score(r, policy) is not None for r in members)
            median = statistics.median(score(r, policy) for r in members) if available else None
            require(available == saved['available_by_policy'][policy] and median == saved['median_R_kcal_mol'][policy], 'Existing triple availability/median differs')
            predictions = calls(median, policy, expected)
            require(predictions['existing_bands'] == saved['decisions'][policy], 'Triple old-band replay differs')
            for rule, decision in predictions.items():
                item = {k: saved[k] for k in ('triple_id', 'protein', 'members', 'pool_case_ids', 'prior_preparation_status')}
                item.update(policy=policy, rule=rule, median_R_kcal_mol=median, expected_class=expected, **decision,
                            old_decision=predictions['existing_bands']['decision'], old_outcome=predictions['existing_bands']['outcome'])
                triples.append(item)
                strata['existing_100_triples', policy, rule].append(item)
    summary = [{'stratum': s, 'policy': p, 'rule': r, 'n': len(v), **counts(v)} for (s, p, r), v in sorted(strata.items())]
    transitions = [r for r in rows if r['rule'] != 'existing_bands' and (r['decision'], r['outcome']) != (r['old_decision'], r['old_outcome'])]
    require(all(t['pass'] for t in tests), 'Real-calibration algebra check failed')
    for name, values in [('CONTEXTS', rows), ('SUMMARY', summary), ('TRANSITIONS', transitions), ('TRIPLES', triples), ('HELD_OUT_CANONICAL', held_out)]:
        table(out / (name + '.tsv'), values)
    result = {'status': 'complete', 'protocol': 'pqq_calibration_only_decision_v1', 'primary_rule': 'canonical_midpoint',
              'models': models, 'summary': summary, 'transitions': transitions, 'tests': tests,
              'held_out_canonical_counts': {p: counts([r for r in held_out if r['policy'] == p]) for p in POLICIES},
              'original_inconclusive_contexts': [r for r in rows if r['old_outcome'] == 'inconclusive'],
              'new_wrong_contexts': [r for r in transitions if r['outcome'] == 'wrong' and r['old_outcome'] != 'wrong'],
              'lost_correct_contexts': [r for r in transitions if r['old_outcome'] == 'correct' and r['outcome'] != 'correct'],
              'denominators': {'unique_contexts': 257, 'memberships': 260, 'canonical_contexts': 25, 'noncanonical_contexts': 232, 'existing_triples': 100, 'superseded_historical_records': 18},
              'historical_ledger': transfer['historical_ledger'], 'superseded_primary_records': transfer['superseded_primary_records'],
              'new_molecular_calls': 0, 'scores_changed': False, 'production_changed': False,
              'limitations': 'Consumed correlated development structures. Sensitivity band is not a confidence interval. Held-out canonical predictions use only their own omitted-group fit; no new blind validation.'}
    write(out / 'RESULT.json', result)
    lines = ['# Calibration-only decision comparison', '', 'Physics and scores unchanged; thresholds use designated canonical25 only.', '', '| Population | Rule | Correct | Wrong | Inconclusive | Unavailable |', '|---|---|---:|---:|---:|---:|']
    for row in summary:
        if row['policy'] == 'operational' and row['stratum'] in ('all_unique', 'canonical', 'noncanonical', 'existing_100_triples'):
            lines.append('| ' + ' | '.join(str(row[k]) for k in ('stratum', 'rule', 'correct', 'wrong', 'inconclusive', 'unavailable')) + ' |')
    lines += ['', 'Every changed call is in TRANSITIONS.tsv. Both selection policies, all source-conditioning strata and per-protein outcomes are in SUMMARY.tsv. Original preparation exclusions and median aggregation remain.', '', result['limitations']]
    (out / 'REPORT.md').write_text('\n'.join(lines) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--inputs-sha', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Numerical comparison requires Slurm')
    placement = subprocess.check_output(['scontrol', 'show', 'job', '-o', os.environ['SLURM_JOB_ID']], text=True)
    fields = dict(v.split('=', 1) for v in placement.split() if '=' in v)
    require(fields.get('JobState') == 'RUNNING' and fields.get('NumCPUs') == '1', 'Exactly one running Slurm CPU required')
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    receipt = {'status': 'failed', 'job_id': os.environ['SLURM_JOB_ID'], 'placement': placement, 'new_molecular_calls': 0}
    try:
        ip = {'path': str(args.inputs.resolve()), 'sha256': args.inputs_sha}
        config = read(ip)
        require(config['implementation'] == pin(__file__) and config['python'] == pin(sys.executable), 'Executable pin changed')
        check(config['plan']); check(config['wrapper'])
        write(out / 'INPUTS.json', {'declaration': ip, **config})
        result = compare(config, out)
        receipt.update(status=result['status'], result=pin(out / 'RESULT.json'))
    except BaseException as error:
        receipt.update(error=str(error), traceback=traceback.format_exc())
    finally:
        receipt['elapsed_seconds'] = time.monotonic() - start
        write(out / 'EXECUTION.json', receipt)
    print(json.dumps(receipt))
    return 0 if receipt['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
