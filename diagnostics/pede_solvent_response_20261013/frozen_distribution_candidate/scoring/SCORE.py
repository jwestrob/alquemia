#!/usr/bin/env python3
"""Saved-result collector; preserves the existing PQQ selection/calibration rules."""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import time

PROTOCOL = 'frozen_distribution_candidate_v1'
ORDER = ('origin', 'adaptive_Ca', 'adaptive_La')
POLICIES = ('operational', 'mathematical', 'origin')
EV = 23.06054783061903
HA = 627.509474


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


@lru_cache(None)
def checked_path(path, sha):
    p = pin(path)
    if p['sha256'] != sha:
        raise ValueError('changed artifact: '+path)
    return json.loads(Path(path).read_text())


def checked(p):
    return checked_path(p['path'], p['sha256'])


def choose(matrix):
    """Same ORDER, strict origin tolerance, and completeness as the pinned baseline.

    The explicit additive solvent field is never mislabeled as a water/vacuum total.
    """
    complete = all(matrix[z].get(q) is not None for z in ('Ca', 'La') for q in ORDER)
    origin_ok = all(matrix[z].get('origin') is not None for z in ('Ca', 'La'))
    selected = {p: None for p in POLICIES}
    for policy in POLICIES:
        if not (origin_ok if policy == 'origin' else complete):
            continue
        mode = {}
        for z in ('Ca', 'La'):
            origin = matrix[z]['origin']
            works = {q: (v['MACE_eV']-origin['MACE_eV'])*EV
                     +(v['solvent_hartree']-origin['solvent_hartree'])*HA
                     for q, v in matrix[z].items() if v is not None}
            best = min(ORDER, key=lambda q: works[q]) if complete else 'origin'
            q = ('origin' if policy == 'origin' else
                 best if policy == 'mathematical' or works[best] < -0.1 else 'origin')
            mode[z] = {'candidate': q, 'work_kcal_mol': works[q],
                       'all_works_kcal_mol': works, 'components': matrix[z][q]}
        ca, la = mode['Ca']['components'], mode['La']['components']
        mode['R_model_kcal_mol'] = ((ca['MACE_eV']-la['MACE_eV'])*EV
                                  +(ca['solvent_hartree']-la['solvent_hartree'])*HA)
        selected[policy] = mode
    return {'status': 'available' if complete else 'unavailable', **selected}


def calibrate(canonical):
    if len(canonical) != 25 or len({r['record_id'] for r in canonical}) != 25:
        raise ValueError('requires exact declared25, including any unavailable members')
    result = {}
    for policy in POLICIES:
        rows = [r for r in canonical if r['pool'][policy] is not None]
        status, gap, bands = 'unavailable_calibration_member', None, None
        if len(rows) == 25:
            maxima = max(r['pool'][policy]['R_model_kcal_mol'] for r in rows if r['expected_class'] == 'Ca')
            minima = min(r['pool'][policy]['R_model_kcal_mol'] for r in rows if r['expected_class'] == 'La')
            gap = minima-maxima
            status = 'available' if gap > 0.02 else 'unsupported_class_separation'
            if status == 'available':
                bands = {'Ca_max': maxima, 'La_min': minima}
        result[policy] = {'status': status, 'available': len(rows), 'gap_kcal_mol': gap, 'bands': bands}
    return result


def decision(score, bands, expected):
    if score is None:
        return {'decision': None, 'outcome': 'unavailable'}
    if bands is None:
        return {'decision': None, 'outcome': 'calibration_unavailable'}
    call = 'Ca' if score <= bands['Ca_max'] else 'La' if score >= bands['La_min'] else 'inconclusive'
    return {'decision': call, 'outcome': 'inconclusive' if call == 'inconclusive'
            else 'correct' if call == expected else 'wrong'}


def decisions(row, calibration):
    return {p: decision(row['pool'][p]['R_model_kcal_mol'] if row['pool'][p] else None,
                        calibration[p]['bands'], row['expected_class']) for p in POLICIES}


def triples(inventory, rows, calibration):
    lookup = {r['case_id']: r for r in rows if r['family'] == 'expanded135'}
    result = []
    for t in inventory['cohorts']['existing_100_triples']:
        members = [lookup.get(k) for k in t['pool_case_ids']]
        available = {p: t['status'] == 'prepared' and all(r is not None and r['pool'][p] is not None
                                                       for r in members) for p in POLICIES}
        medians = {p: statistics.median(r['pool'][p]['R_model_kcal_mol'] for r in members)
                   if available[p] else None for p in POLICIES}
        expected = next((r['expected_class'] for r in members if r is not None), None)
        result.append({'triple_id': t['triple_id'], 'protein': t['protein_id'],
                       'members': t['members'], 'pool_case_ids': t['pool_case_ids'],
                       'prior_preparation_status': t['status'], 'available_by_policy': available,
                       'status': 'available' if available['operational'] else 'unavailable',
                       'median_R_kcal_mol': medians,
                       'decisions': {p: decision(medians[p], calibration[p]['bands'], expected) for p in POLICIES}})
    return result


def score_row(record, native_tasks, baseline_calibration):
    matrix = {z: {} for z in ('Ca', 'La')}
    failures = []
    for z in matrix:
        for q in ORDER:
            cell = record['matrix'][z][q]
            task_id = cell.get('native_task_id')
            native = native_tasks.get(task_id)
            reason = cell.get('reason') if task_id is None else 'native_cell_not_collected'
            matrix[z][q] = None
            if native is not None:
                reason = native.get('reason')
                if native['status'] == 'available':
                    checked(native['result'])
                    g = native['G_qvac_hartree']
                    mace = checked(cell['MACE'])
                    if not math.isfinite(g) or mace['energy_eV'] != cell['MACE_eV']:
                        raise ValueError('invalid accepted native/MACE energy')
                    matrix[z][q] = {'MACE_eV': cell['MACE_eV'], 'solvent_hartree': g,
                                    'solvent_definition': 'native_ALPB_at_GFN2_vacuum_charges',
                                    'native_task_id': task_id, 'native_result': native['result'],
                                    'xyz': cell['xyz'], 'MACE': cell['MACE']}
            if matrix[z][q] is None:
                failures.append({'metal': z, 'candidate': q, 'task_id': task_id, 'reason': reason})
    old_collection = checked(record['baseline_collection'])
    old = old_collection['rows'][record['baseline_collection_row_index']]
    if old['case_id'] != record['case_id'] or old['source_case_id'] != record['source_case_id']:
        raise ValueError('baseline source join differs')
    row = {k: record[k] for k in ('record_id', 'case_id', 'source_case_id', 'selection_id', 'protein',
                                  'biological_group', 'role', 'expected_class', 'source_conditioning_metal', 'family')}
    row.update(pool=choose(matrix), missing_cells=failures,
               baseline={'collection': record['baseline_collection'], 'row_index': record['baseline_collection_row_index'],
                         'pool': old['pool']})
    row['baseline']['decisions'] = decisions({**row, 'pool': old['pool']}, baseline_calibration)
    return row


def summarize(ids, lookup):
    rows = [lookup[k] for k in ids]
    result = {'declared': len(ids), 'policies': {}}
    for p in POLICIES:
        both = [r for r in rows if r['pool'][p] is not None and r['baseline']['pool'][p] is not None]
        result['policies'][p] = {'candidate': dict(Counter(r['decisions'][p]['outcome'] for r in rows)),
            'baseline': dict(Counter(r['baseline']['decisions'][p]['outcome'] for r in rows)),
            'common_coverage': len(both),
            'common_candidate': dict(Counter(r['decisions'][p]['outcome'] for r in both)),
            'common_baseline': dict(Counter(r['baseline']['decisions'][p]['outcome'] for r in both)),
            'transitions': dict(Counter(r['baseline']['decisions'][p]['outcome']+'->'+r['decisions'][p]['outcome'] for r in rows))}
    return result


def run(config, stage, native_collection, canonical_calibration=None):
    inventory = checked(config['inventory'])
    canonical_ids = inventory['calibration_rule']['designated_record_ids']
    if config['protocol'] != PROTOCOL or len(canonical_ids) != 25:
        raise ValueError('candidate protocol or declared canonical population differs')
    native = checked(native_collection)
    if (native['status'] != 'complete' or native['stage'] != stage
            or native['manifest'] != config['native_manifest'] or native['inventory'] != config['inventory']):
        raise ValueError('native collection identity/stage/status differs')
    native_tasks = dict(native['tasks'])
    if stage == 'transfer':
        calibration = checked(canonical_calibration)
        if (calibration['protocol'] != PROTOCOL or calibration['inventory'] != config['inventory']
                or calibration['native_manifest'] != config['native_manifest']
                or calibration['calibrations']['operational']['status'] != 'available'
                or calibration['scorer'] != config['scorer']):
            raise ValueError('transfer requires exact frozen passed canonical calibration')
        primary_native = checked(calibration['canonical_native_collection'])
        for k, v in primary_native['tasks'].items():
            if k in native_tasks and native_tasks[k] != v:
                raise ValueError('canonical native result silently changed during transfer')
            native_tasks[k] = v
        calibrations = calibration['calibrations']
    baseline_calibration = checked(inventory['baseline_canonical'])['calibration']
    records = [r for r in inventory['records'] if stage == 'transfer' or r['record_id'] in canonical_ids]
    rows = [score_row(r, native_tasks, baseline_calibration) for r in records]
    lookup = {r['record_id']: r for r in rows}
    if stage == 'canonical':
        calibrations = calibrate([lookup[k] for k in canonical_ids])
    else:
        for old_row in calibration['rows']:
            if lookup[old_row['record_id']]['pool'] != old_row['pool']:
                raise ValueError('canonical scores changed after bands frozen')
    for row in rows:
        row['decisions'] = decisions(row, calibrations)
        row['old_band_transfer'] = decisions(row, baseline_calibration)
    groups = {'canonical25': summarize(canonical_ids, lookup)}
    grouped_triples = []
    if stage == 'transfer':
        for name in ('expanded135', 'primary_ca125', 'two_ca15'):
            groups[name] = summarize(inventory['cohorts'][name], lookup)
        overlay = inventory['cohorts']['current_ca125_annotated_overlay']
        groups['current_ca125_annotated_overlay'] = summarize([r['record_id'] for r in overlay], lookup)
        grouped_triples = triples(inventory, rows, calibrations)
    return {'protocol': PROTOCOL, 'stage': stage, 'status': 'complete', 'inventory': config['inventory'],
            'native_manifest': config['native_manifest'], 'scorer': config['scorer'],
            'canonical_native_collection': native_collection if stage == 'canonical' else calibration['canonical_native_collection'],
            'canonical_calibration_receipt': canonical_calibration, 'native_collection': native_collection,
            'calibrations': calibrations, 'baseline_calibrations': baseline_calibration,
            'rows': rows, 'cohorts': groups, 'triples': grouped_triples,
            'current_ca125_annotated_overlay': inventory['cohorts']['current_ca125_annotated_overlay'],
            'production_changed': False, 'new_MACE_SCC_DFT_calls': 0,
            'interpretation': 'consumed development references; separately calibrated fixed-distribution candidate; no affinity claim'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True)
    p.add_argument('--sha', required=True)
    p.add_argument('--stage', choices=('canonical', 'transfer'), required=True)
    p.add_argument('--native-collection', required=True)
    p.add_argument('--canonical-calibration')
    p.add_argument('--output', required=True)
    a = p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Slurm allocation required')
    config = checked({'path': str(Path(a.config).resolve()), 'sha256': a.sha})
    if pin(__file__) != config['scorer']:
        raise SystemExit('unfrozen scorer implementation')
    started = time.monotonic()
    result = run(config, a.stage, pin(a.native_collection), pin(a.canonical_calibration) if a.canonical_calibration else None)
    result.update(elapsed_seconds=time.monotonic()-started, job_id=os.environ['SLURM_JOB_ID'], config=pin(a.config))
    with Path(a.output).open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({'stage': a.stage, 'calibrations': result['calibrations'], 'cohorts': result['cohorts']}, indent=2))


if __name__ == '__main__':
    main()
