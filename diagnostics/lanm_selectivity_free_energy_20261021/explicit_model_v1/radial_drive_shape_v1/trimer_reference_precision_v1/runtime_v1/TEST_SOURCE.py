#!/usr/bin/env python3
"""Focused saved/ manufactured controls. Every real Popen is forbidden."""
import copy
import json
from pathlib import Path
import shutil
import time
from unittest.mock import patch
import RUN as run

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p/'diagnostics').is_dir())
WORK = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/radial_drive_shape_v1/trimer_reference_precision_v1/runtime_v1'
R = run.residuals


def manufactured(algorithm='SOSCF', energy='-1.0000e-12', rotation='1.0000e-10', density='2.0000e-12'):
    diag = '2.0000e-03' if algorithm == 'SOSCF' else '3.0000e-10'
    rows = ('----------------S-O-S-C-F----------------\n'
            '1 -100.0 1.00e-12 1.00e-12 '+density+' 3.00e-10 1.0\n') if algorithm == 'SOSCF' else (
            '----------------D-I-I-S----------------\n'
            '1 -100.0 1.00e-12 1.00e-12 '+density+' 3.00e-10 0.0 1.0\n')
    return ('Convergence Check Mode ConvCheckMode .... All-Criteria\n'
            'Convergence forced ConvForced .... 1\nMaximum # iterations MaxIter .... 300\n'
            'Energy Change TolE .... 1.000e-11 Eh\nOrbital Gradient TolG .... 1.000e-08\n'
            'Orbital Rotation angle TolX .... 1.000e-08\nDIIS Error TolErr .... 1.000e-08\n'
            'ORCA LEAN-SCF\n'+rows+'SCF CONVERGED AFTER 1 CYCLES\nSCF CONVERGENCE\n'
            'Last Energy change ... '+energy+' Tolerance : 1.0000e-11\n'
            'Last MAX-Density change ... '+density+' Tolerance : 1.0000e-08\n'
            'Last RMS-Density change ... 1.0000e-12 Tolerance : 1.0000e-09\n'
            'Last DIIS Error ... '+diag+' Tolerance : 1.0000e-08\n'
            'Last Orbital Gradient ... 3.0000e-10 Tolerance : 1.0000e-08\n'
            'Last Orbital Rotation ... '+rotation+' Tolerance : 1.0000e-08\nORBITAL ENERGIES\n')


def main():
    tick, cpu = time.monotonic(), time.process_time()
    out = WORK/'controls_v1'
    out.mkdir(parents=True, exist_ok=False)
    checks = []
    result = {'status': 'started', 'implementation': run.pin(__file__), 'native_calls': 0,
              'real_subprocesses': 0, 'error': None, 'manufactured_results_are_not_reference_admissions': True}

    def check(label, ok):
        run.require(ok, label)
        checks.append(label)

    def reject(label, f):
        try:
            f()
        except (ValueError, FileExistsError):
            checks.append(label)
        else:
            raise AssertionError('Did not reject: '+label)

    try:
        with patch.object(run.subprocess, 'Popen', side_effect=AssertionError('Real process forbidden')):
            cfg = json.loads((WORK/'preparation_v2/CONFIG_SCOUT_TEMPLATE.json').read_text())
            m, inv, d, common, collector, qhelper = run.context(cfg)
            result['manifest'] = cfg['manifest']
            check('Exact prepared source/decks/helper context passes without a process', len(m['rows']) == 14)
            for row in m['rows']:
                old = run.verify(row['old']['input']).read_text()
                check('SCF-only inverse: '+row['cell_id'], run.verify(row['input']).read_text().replace(run.NEW_SCF, run.OLD_SCF) == old)
            reject('Repeated SCF block rejected', lambda: run.deck(run.OLD_SCF+run.OLD_SCF))
            old_row = next(r for r in m['rows'] if r['cell_id'] == d['staging']['scout_order'][0])
            old_raw = run.verify(old_row['old']['output']).read_text()
            reject('Actual old mode2 output cannot pass new protocol', lambda: R.current(old_raw))
            check('Manufactured SOSCF accepts stale DIIS only descriptively', R.current(manufactured())['all_required_passed'])
            check('Manufactured opposite energy-summary sign passes magnitude join', R.current(manufactured(energy='1.0000e-12'))['all_required_passed'])
            check('Manufactured terminal DIIS joins its current error', R.current(manufactured('DIIS'))['all_required_passed'])
            check('Current rotation gate rejects large value', not R.current(manufactured(rotation='1.0000e-04'))['all_required_passed'])
            check('Current density gate rejects matching-but-large value', not R.current(manufactured(density='2.0000e-05'))['all_required_passed'])
            bad_energy = manufactured(energy='-2.0000e-08').replace('1.00e-12 1.00e-12', '2.00e-08 1.00e-12', 1)
            check('Negative large deltaE never passes', not R.current(bad_energy)['all_required_passed'])
            for label, text in [
                ('Wrong mode', manufactured().replace('All-Criteria', 'Total+1el-Energy')),
                ('Wrong forced convergence', manufactured().replace('ConvForced .... 1', 'ConvForced .... 0')),
                ('Wrong declared threshold', manufactured().replace('TolE .... 1.000e-11', 'TolE .... 1.000e-09')),
                ('Missing rotation summary', manufactured().replace('Last Orbital Rotation', 'Missing Orbital Rotation')),
                ('Mismatched final gradient', manufactured().replace('Last Orbital Gradient ... 3.0000e-10', 'Last Orbital Gradient ... 9.0000e-10')),
                ('Mismatched final count', manufactured().replace('AFTER 1 CYCLES', 'AFTER 2 CYCLES')),
                ('Unknown solver', manufactured().replace('S-O-S-C-F', 'T-R-A-H')),
                ('Duplicate summary', manufactured()+'SCF CONVERGENCE\n'),
                ('Nonfinite residual', manufactured().replace('Last Orbital Rotation ... 1.0000e-10', 'Last Orbital Rotation ... NaN')),
                ('Negative density', manufactured(density='-2.0000e-12'))]:
                reject(label+' rejected', lambda text=text: R.current(text))
            check('Printing join uses literal precision', R.print_join('3.8286e-7', '3.83e-7') and not R.print_join('3.8286e-7', '3.80e-7'))
            check('Old settings signature is exact', R.settings(old_raw) == old_row['baseline_grid_integral_settings'])
            check('Grid change changes signature', R.settings(old_raw.replace('1.0000e-12', '2.0000e-12')) != old_row['baseline_grid_integral_settings'])
            check('Printed occupations retain actual electron total', sum(R.dec(v) for _,v in R.occupations(old_raw)) == old_row['old']['cell']['expected_electrons'])
            # Saved original native identity parses unchanged; it remains ineligible for the new residual protocol.
            parsed = common.helper.parse(old_raw, old_row['old']['cell'], 0)
            check('Saved original PBE0 identity still passes', parsed['status'] == 'complete' and all(parsed['checks'].values()))
            # Pure release tests use in-memory synthetic objects, never a real execution release file.
            cp = {'path': 'manufactured_config', 'sha256': 'synthetic'}
            rp = {'path': 'manufactured_release', 'sha256': 'synthetic'}
            sp = {'path': 'manufactured_review', 'sha256': 'synthetic'}
            rel = json.loads((WORK/'preparation_v2/ROOT_RELEASE_TEMPLATE.json').read_text())
            rel.update(configuration=cp, source_review=sp, login_reservation_confirmed=True,
                       other_owned_login_numerical_threads_at_release=0)
            review = {'status': 'passed_reference_precision_runtime_source_review', 'implementation': run.pin(HERE/'RUN.py'),
                      'manifest': cfg['manifest'], 'residual_parser': m['residual_parser']}
            def fake_read(p):
                return rel if p == rp else review if p == sp else (_ for _ in ()).throw(ValueError('Unexpected synthetic pin'))
            with patch.object(run, 'read', fake_read):
                check('Exact two-state scout release accepted by pure contract', run.release(cfg, cp, rp, m, d)[1] == d['staging']['scout_order'])
                rel['cell_ids'] = rel['cell_ids']+['extra']
                reject('Extra state release rejected', lambda: run.release(cfg, cp, rp, m, d))
                rel['cell_ids'] = d['staging']['scout_order']
                rel['other_owned_login_numerical_threads_at_release'] = 9
                reject('Overcommitted login reservation rejected', lambda: run.release(cfg, cp, rp, m, d))
            # Two tiny artificial state machines exercise durable intent and fail-stop with no physical readback.
            for fail_first in (False, True):
                base = out/('manufactured_failure' if fail_first else 'manufactured_success')
                fake_cfg = dict(cfg, campaign_root=str(base/'campaign'), output_root=str(base/'campaign/scout'))
                base.mkdir()
                run.save(base/'config.json', fake_cfg)
                run.save(base/'release.json', {'synthetic_nonexecuting_fixture': True})
                fake_m = dict(m, campaign_root=fake_cfg['campaign_root'], rows=[r for r in m['rows'] if r['cell_id'] in d['staging']['scout_order']])
                calls = []
                def fake_invoke(argv, folder, env):
                    check('Attempt intent exists before mocked invocation', (folder/'START.json').exists())
                    calls.append(folder.name)
                    receipt = {'returncode': 0, 'synthetic': True}
                    run.save(folder/'NATIVE_RECEIPT.json', receipt)
                    return receipt
                def fake_admit(*args):
                    if fail_first:
                        raise ValueError('Manufactured residual failure')
                    return {'synthetic': True}, {'synthetic': True}, {'synthetic': True}
                with patch.object(run, 'context', return_value=(fake_m, inv, d, None, None, None)), \
                     patch.object(run, 'release', return_value=({}, d['staging']['scout_order'], None)), \
                     patch.object(run, 'runtime', return_value={'synthetic': True}), \
                     patch.object(run, 'invoke', fake_invoke), patch.object(run, 'admit', fake_admit):
                    rr = run.run(base/'config.json', base/'release.json')
                    check('Failure stops after one attempt' if fail_first else 'Exactly two serial mocked attempts', len(calls) == (1 if fail_first else 2))
                    check('Failure retained' if fail_first else 'Scout terminal denominator retained',
                          rr['status'] == ('reference_precision_phase_failed' if fail_first else 'reference_precision_scout_complete'))
                    reject('Existing phase output cannot execute again', lambda: run.run(base/'config.json', base/'release.json'))
            result['status'] = 'passed_reference_precision_runtime_controls'
    except Exception as exc:
        result.update(status='failed_reference_precision_runtime_controls', error=repr(exc))
    result.update(checks=checks, check_count=len(checks), runtime_implementation=run.pin(HERE/'RUN.py'),
                  residual_parser=run.pin(HERE/'RESIDUALS.py'),
                  cost={'wall_seconds': time.monotonic()-tick, 'CPU_seconds': time.process_time()-cpu, 'numerical_threads': 1})
    run.save(out/'RESULT.json', result)
    print(json.dumps({k: result[k] for k in ('status', 'error', 'check_count', 'cost')}))
    return 0 if result['error'] is None else 1


if __name__ == '__main__':
    raise SystemExit(main())
