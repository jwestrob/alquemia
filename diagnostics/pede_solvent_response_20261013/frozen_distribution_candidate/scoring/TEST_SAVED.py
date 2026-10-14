#!/usr/bin/env python3
"""Regression of unchanged policy on actual archived baseline energies, not new chemistry.

Missing-input checks deliberately corrupt in-memory copies of those real fixtures.
"""
import argparse
import copy
import json
import math
import os
from pathlib import Path
import time
import SCORE as s


def equal_pool(a, b):
    assert a['status'] == b['status']
    for p in s.POLICIES:
        assert (a[p] is None) == (b[p] is None)
        if a[p] is None:
            continue
        assert math.isclose(a[p]['R_model_kcal_mol'], b[p]['R_model_kcal_mol'], abs_tol=1e-7, rel_tol=0)
        for z in ('Ca', 'La'):
            assert a[p][z]['candidate'] == b[p][z]['candidate']
            assert math.isclose(a[p][z]['work_kcal_mol'], b[p][z]['work_kcal_mol'], abs_tol=1e-7, rel_tol=0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--inventory', required=True)
    p.add_argument('--sha', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Slurm required for real-fixture analysis')
    started = time.monotonic()
    inventory_pin = {'path': str(Path(a.inventory).resolve()), 'sha256': a.sha}
    inv = s.checked(inventory_pin)
    rule = inv['calibration_rule']
    assert rule['requires_all'] == 25 and rule['minimum_gap_strictly_greater_kcal_mol'] == 0.02
    assert inv['selection_rule']['candidate_order'] == list(s.ORDER)
    assert inv['selection_rule']['origin_retained_unless_minimum_work_below_kcal_mol'] == -0.1
    rows, matrices, original = [], {}, {}
    for r in inv['records']:
        old = s.checked(r['baseline_collection'])['rows'][r['baseline_collection_row_index']]
        assert old['case_id'] == r['case_id'] and old['source_case_id'] == r['source_case_id']
        matrix = {z: {} for z in ('Ca', 'La')}
        for z in matrix:
            for q in s.ORDER:
                c = old.get('components', {}).get(z+'__'+q)
                cell = r['matrix'][z][q]
                matrix[z][q] = None
                if c and c.get('accepted') and cell['status'] == 'available_archived_inputs':
                    # Real archived COMPLETE water-minus-vacuum contribution, used
                    # solely to regression-test the unchanged algebra/policy.
                    matrix[z][q] = {'MACE_eV': cell['MACE_eV'],
                        'solvent_hartree': c['primary_water_energy_hartree']-c['primary_vacuum_energy_hartree']}
        row = {**r, 'pool': s.choose(matrix)}
        equal_pool(row['pool'], old['pool'])
        rows.append(row)
        matrices[r['record_id']] = matrix
        original[r['record_id']] = old
    lookup = {r['record_id']: r for r in rows}
    canonical = [lookup[k] for k in rule['designated_record_ids']]
    calibration = s.calibrate(canonical)
    baseline = s.checked(inv['baseline_canonical'])['calibration']
    for p in s.POLICIES:
        assert calibration[p]['status'] == baseline[p]['status'] == 'available'
        for k in ('Ca_max', 'La_min'):
            assert math.isclose(calibration[p]['bands'][k], baseline[p]['bands'][k], rel_tol=0, abs_tol=1e-7)
        assert s.decision(baseline[p]['bands']['Ca_max'], baseline[p]['bands'], 'Ca')['decision'] == 'Ca'
        assert s.decision(baseline[p]['bands']['La_min'], baseline[p]['bands'], 'La')['decision'] == 'La'
    for row in rows:
        calculated = s.decisions(row, baseline)
        for p in s.POLICIES:
            if 'decisions' in original[row['record_id']]:
                assert calculated[p] == original[row['record_id']]['decisions'][p]
    grouped = s.triples(inv, rows, baseline)
    old_grouped = s.checked(inv['sources']['expanded135']['collection'])['triples']
    assert len(grouped) == len(old_grouped) == 100
    for t, old in zip(grouped, old_grouped):
        assert t['triple_id'] == old['triple_id'] and t['pool_case_ids'] == old['pool_case_ids']
        assert t['available_by_policy'] == old['available_by_policy'] and t['decisions'] == old['decisions']
        for p in s.POLICIES:
            x, y = t['median_R_kcal_mol'][p], old['median_R_kcal_mol'][p]
            assert (x is None) == (y is None)
            if x is not None:
                assert math.isclose(x, y, rel_tol=0, abs_tol=1e-7)
    # Explicitly corrupted REAL fixtures: never presented as molecular evidence.
    broken = copy.deepcopy(matrices[canonical[0]['record_id']])
    broken['La']['adaptive_Ca'] = None
    failed = s.choose(broken)
    assert failed['status'] == 'unavailable' and failed['operational'] is None and failed['mathematical'] is None
    assert failed['origin'] is not None
    missing_cal = copy.deepcopy(canonical)
    missing_cal[0]['pool'] = failed
    missing = s.calibrate(missing_cal)
    assert missing['operational']['status'] == 'unavailable_calibration_member' and missing['operational']['bands'] is None
    assert s.decision(None, baseline['operational']['bands'], 'Ca')['outcome'] == 'unavailable'
    assert s.decision(canonical[0]['pool']['origin']['R_model_kcal_mol'], None, 'Ca')['outcome'] == 'calibration_unavailable'
    report = {'status': 'passed', 'inventory': inventory_pin, 'scorer': s.pin(s.__file__),
              'tests_implementation': s.pin(__file__), 'job_id': os.environ['SLURM_JOB_ID'],
              'archived_pools_checked': len(rows), 'canonical_members_checked': len(canonical),
              'archived_triples_checked': len(grouped), 'missing_input_checks': 'deliberately corrupted copies of real fixtures',
              'new_molecular_calls': 0, 'candidate_scores_calculated': 0,
              'elapsed_seconds': time.monotonic()-started}
    with Path(a.output).open('x') as f:
        json.dump(report, f, indent=2)
        f.write('\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
