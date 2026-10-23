"""Read the 68 completed water states once; apply unchanged 32 FD checks.

Saved-data analysis only, on Slurm standard_expand. Reuses frozen COLLECT.slope;
does not launch native calls, relax tolerances, or admit the donor regression.
"""
import argparse
import json
import os
import time
from pathlib import Path
from PREPARE import pin, read, save
from COLLECT import slope

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scout', type=Path, required=True)
    ap.add_argument('--partial', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    t, cpu = time.monotonic(), time.process_time()
    scout, partial = read(pin(a.scout)), read(pin(a.partial))
    assert scout['status'] == 'native_force_scout_complete'
    assert partial['status'] == 'failed'
    assert partial['reused_scout_rows'] == scout['rows']
    panel = read(scout['panel'])
    rows = {}
    for state, authority in {**scout['rows'], **partial['rows']}.items():
        if state.startswith('aquo'):
            row = read(authority)
            assert row['status'] == 'native_force_state_complete'
            assert row['state'] == panel['states'][state]
            rows[state] = row
    assert len(rows) == 68
    checks = {}
    for state, row in rows.items():
        s = row['state']
        if s['sign'] != 1:
            continue
        center = rows[s['case_id'] + '__' + s['coefficient_mode'] + '__center']
        minus = rows[state[:-4] + 'minus']
        checks[state] = slope(row['energy_kcal_mol'], minus['energy_kcal_mol'],
                              s['step_A'], center['gradient_kcal_mol_A'],
                              panel['directions'][s['case_id']][s['direction']])
    assert len(checks) == 32
    passed = all(c['passed'] for c in checks.values())
    result = dict(status='passed_saved_water_force_checks' if passed else 'failed_directional_gate',
                  implementation=pin(Path(__file__).resolve()),
                  collector=pin(Path(__file__).with_name('COLLECT.py').resolve()),
                  scout=pin(a.scout.resolve()), partial=pin(a.partial.resolve()),
                  panel=scout['panel'], checks=checks, reused_states=68,
                  new_native_calls=0, original_donor_regression_unresolved=True,
                  job_id=os.environ['SLURM_JOB_ID'],
                  cost=dict(wall_seconds=time.monotonic()-t, CPU_seconds=time.process_time()-cpu))
    a.out.mkdir(parents=True, exist_ok=False)
    save(a.out/'RESULT.json', result)
    print(json.dumps({'status':result['status'], 'checks':len(checks)}))
    return 0 if passed else 1

if __name__ == '__main__':
    raise SystemExit(main())
