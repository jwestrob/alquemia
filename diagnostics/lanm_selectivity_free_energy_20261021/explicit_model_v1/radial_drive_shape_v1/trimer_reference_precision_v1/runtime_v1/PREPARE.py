#!/usr/bin/env python3
"""Freeze14 SCF-only input changes and source/runtime authorities; zero native."""
import json
import os
from pathlib import Path
import time
import RUN as run

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p/'diagnostics').is_dir())
WORK = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/radial_drive_shape_v1/trimer_reference_precision_v1/runtime_v1'


def main():
    tick, cpu = time.monotonic(), time.process_time()
    out = WORK/'preparation_v1'
    out.mkdir(parents=True, exist_ok=False)
    result = {'status': 'started', 'implementation': run.pin(__file__), 'native_calls': 0, 'error': None}
    try:
        x = run.read(run.pin(HERE/'METHOD_INPUTS.json'))
        d = run.read(x['declaration'])
        root = run.read(x['protocol_review'])
        run.require(root['status'] == 'admitted_reference_precision_protocol_design' and root['declaration'] == x['declaration'], 'Root protocol differs')
        run.read(root['independent_review'])
        inv = run.read(d['inventory'])
        expected = {'ConvCheckMode': 0, 'ConvForced': 1, 'TolE': '1e-11', 'TolMAXP': '1e-8',
                    'TolRMSP': '1e-9', 'TolErr': '1e-8', 'TolG': '1e-8', 'TolX': '1e-8'}
        run.require(d['frozen_scf_termination'] == expected and d['new_reference_inputs_allowed_for_preparation'] == 14,
                    'Exact preparation authority differs')
        for p in x['methods'].values():
            run.verify(p)
        inputs = out/'inputs'
        inputs.mkdir()
        rows = []
        for old in inv['rows']:
            original = run.verify(old['input']).read_text()
            new = run.deck(original)
            run.require(new.replace(run.NEW_SCF, run.OLD_SCF) == original, 'SCF-only inverse mismatch')
            target = inputs/(old['cell_id']+'.inp')
            target.write_text(new)
            raw = run.verify(old['output']).read_text()
            occ = run.residuals.occupations(raw)
            run.require(len(occ) == old['cell']['expected_basis_functions'], 'Original occupation count differs')
            rows.append({'cell_id': old['cell_id'], 'old': old, 'input': run.pin(target),
                         'baseline_grid_integral_settings': run.residuals.settings(raw), 'baseline_occupations': occ,
                         'only_changed_block': '%scf', 'inverse_matches_exact_original': True})
        m = {'status': 'prepared_reference_precision_login_inputs', 'preparer': run.pin(__file__),
             'implementation': run.pin(HERE/'RUN.py'), 'residual_parser': run.pin(HERE/'RESIDUALS.py'),
             'method_inputs': run.pin(HERE/'METHOD_INPUTS.json'), 'declaration': x['declaration'],
             'protocol_review': x['protocol_review'], 'inventory': d['inventory'], 'methods': x['methods'],
             'orca': inv['orca'], 'runtime_environment': x['runtime_environment'], 'rows': rows,
             'maximum_attempts': 14, 'scout_ids': d['staging']['scout_order'],
             'remaining_ids': d['staging']['continuation_order'], 'campaign_root': str(WORK/'run_v1/execution'),
             'new_native_calls': 0, 'execution_release': None}
        run.save(out/'MANIFEST.json', m)
        cfg = {'status': 'bound_reference_precision_login_configuration', 'implementation': m['implementation'],
               'manifest': run.pin(out/'MANIFEST.json'), 'phase': 'scout', 'campaign_root': m['campaign_root'],
               'output_root': str(Path(m['campaign_root'])/'scout'), 'scout_result': None, 'scout_root_review': None}
        run.save(out/'CONFIG_SCOUT_TEMPLATE.json', cfg)
        remaining = dict(cfg, phase='remaining', output_root=str(Path(m['campaign_root'])/'remaining'))
        run.save(out/'CONFIG_REMAINING_TEMPLATE.json', remaining)
        release = {'status': 'approved_reference_precision_login_phase', 'configuration': None,
                   'implementation': m['implementation'], 'manifest': cfg['manifest'], 'declaration': m['declaration'],
                   'protocol_review': m['protocol_review'], 'source_review': None, 'phase': 'scout',
                   'cell_ids': m['scout_ids'], 'maximum_phase_attempts': 2, 'maximum_total_attempts': 14,
                   'campaign_root': m['campaign_root'], 'execution_host': 'biotite', 'scheduler_job_id': None,
                   'reserved_login_numerical_threads': 8, 'total_owned_login_limit': 16,
                   'login_reservation_confirmed': False, 'other_owned_login_numerical_threads_at_release': None,
                   'scout_result': None, 'scout_root_review': None}
        run.save(out/'ROOT_RELEASE_TEMPLATE.json', release)
        result.update(status='prepared_reference_precision_runtime_source', manifest=cfg['manifest'],
                      inputs_prepared=len(rows), exact_inverse_checks=len(rows), original_occupation_counts_checked=len(rows),
                      original_grid_integral_signatures_copied=len(rows), controls_and_source_review_pending=True)
    except Exception as exc:
        result.update(status='failed_reference_precision_runtime_preparation', error=repr(exc))
    result['cost'] = {'wall_seconds': time.monotonic()-tick, 'CPU_seconds': time.process_time()-cpu,
                      'numerical_threads': 1, 'PID': os.getpid(), 'scope': 'Metadata/source/deck preparation only'}
    run.save(out/'RESULT.json', result)
    print(json.dumps({k: result[k] for k in ('status', 'error', 'cost')}))
    return 0 if result['error'] is None else 1


if __name__ == '__main__':
    raise SystemExit(main())
