#!/usr/bin/env python3
"""Archive terminal metadata for the bounded failed pilot; no scientific calls."""
import csv
import datetime
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
D = Path(__file__).resolve().parent
W = ROOT / 'workspaces/lady_electronic_decomposition_20261023'
TASK = W / 'cp2k_quartet_memory_v1/tasks/aquo8_central_La'
CONFIG = D / 'cp2k_memory_recovery_v1/CONFIG.json'
LEDGER = D / 'cp2k_inputs_v1/ENERGY_LEDGER.json'
GOAL = Path('/home/jwestrob/.codex/attachments/2ef9917f-aee5-4f5d-ba3a-af8ff24784c7/pasted-text-1.txt')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def main():
    assert not (D / 'FINAL_RESULT.json').exists()
    cfg = json.loads(CONFIG.read_text())
    ledger = json.loads(LEDGER.read_text())
    run = json.loads((TASK / 'RUN_RESULT.json').read_text())
    native = run['commands'][1]
    assert run['job_id'] == '21623' and run['status'] == 'failed'
    assert native['timed_out'] is True and native['cap_seconds'] == 7073
    assert run['config_sha256'] == sha(CONFIG)
    pins = []
    for pin in cfg['pins']:
        actual = sha(pin['path'])
        assert actual == pin['sha256'], pin['path']
        pins.append({'path': pin['path'], 'sha256': actual})
    assert not Path(cfg['pilot_admission']).exists()
    assert not (TASK / 'RESULT.json').exists()  # Timeout precedes parser invocation.
    assert len(run['commands']) == 2
    rows = (TASK / 'execution/endpoint.out').read_text().splitlines()
    ot = [{'line': i + 1, 'raw': line} for i, line in enumerate(rows)
          if re.match(r'^\s+\d+ OT CG\s+', line)]
    assert len(ot) == 43
    markers = {label: sum(token in line for line in rows) for label, token in {
        'fragment_converged': 'SCF run converged',
        'R_fragment_sum': 'Single-molecule energy:',
        'F_frozen': 'Energy of the initial guess:',
        'P_block_trace': 'ALMO SCF DIIS',
        'U_full_energy': 'Total FORCE_EVAL',
        'normal_end': 'PROGRAM ENDED AT',
    }.items()}
    assert not any(markers.values()), markers
    mirrored = []
    dest = Path(cfg['diagnostic_results']) / cfg['pilot']
    for artifact in run['artifacts']:
        src = Path(artifact['path'])
        mirror = dest / src.name
        assert sha(src) == artifact['sha256'] == sha(mirror)
        mirrored.append({'path': str(mirror.relative_to(D)), 'sha256': sha(mirror)})
    assert sha(TASK / 'RUN_RESULT.json') == sha(dest / 'RUN_RESULT.json')
    for name in ['SUBMISSION.json', 'completion_event.json', 'job.sbatch']:
        target = dest / name
        assert not target.exists()
        shutil.copyfile(TASK / name, target)
        mirrored.append({'path': str(target.relative_to(D)), 'sha256': sha(target)})
    shutil.copyfile(GOAL, D / 'GOAL_SCOPE.txt')
    command = ['sacct', '-j', '21396,21440,21618,21623',
               '--format=JobID,State,ExitCode,ElapsedRaw,AllocCPUS,ReqMem,MaxRSS', '-P']
    account = subprocess.run(command, check=True, text=True, capture_output=True).stdout
    records = list(csv.DictReader(io.StringIO(account), delimiter='|'))
    parents = [r for r in records if '.' not in r['JobID']]
    assert len(parents) == 4 and all(r['State'] in ['FAILED', 'COMPLETED', 'OUT_OF_MEMORY'] for r in parents)
    assert next(r for r in parents if r['JobID'] == '21623')['State'] == 'FAILED'
    allocated = sum(int(r['ElapsedRaw']) * int(r['AllocCPUS']) for r in parents)
    prior = json.loads((W / 'cp2k_quartet_v1/tasks/aquo8_central_La/RUN_RESULT.json').read_text())
    stamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    costs = {'observed_utc': stamp, 'source_command': command, 'scheduler_rows': records,
             'terminal_allocated_cpu_seconds': allocated, 'running_job_cost_excluded_from_terminal_sum': [],
             'kind': 'allocated_capacity_not_measured_utilization', 'new_admitted_scientific_cells': 0,
             'setup_job_ids': ['21396', '21440'], 'pilot_attempt_job_ids': ['21618', '21623'],
             'native_cap_seconds_recovery': 7073,
             'native_elapsed_seconds_including_termination': native['elapsed_seconds'],
             'total_native_elapsed_seconds_two_attempts': prior['commands'][1]['elapsed_seconds'] + native['elapsed_seconds'],
             'no_new_submissions_by_finalization': True}
    shutil.copyfile(D / 'COSTS.json', D / 'COSTS_BEFORE_TERMINAL.json')
    (D / 'COSTS.json').write_text(json.dumps(costs, indent=2) + '\n')
    scope = []
    for gid, system in cfg['systems'].items():
        task = Path(system['task'])
        if gid != cfg['pilot']:
            assert not (task / 'SUBMISSION.json').exists() and not (task / 'RUN_RESULT.json').exists()
        assert sha(system['xyz']) == system['xyz_sha256']
        assert sha(system['deck']) == system['deck_sha256']
        scope.append({'geometry_id': gid, 'status': 'pilot_timeout' if gid == cfg['pilot'] else 'not_run_pilot_not_admitted',
                      'xyz_sha256': system['xyz_sha256'], 'input_sha256': system['deck_sha256']})
    dd = ledger['double_difference']
    table = [
        ('archived_full_PBE0_D4', str(dd['published_reference_value_kcal_mol']), 'kcal/mol', 'archived_context'),
        ('archived_electronic_PBE0', dd['electronic_kcal_mol_from_tokens'], 'kcal/mol', 'archived_context'),
        ('archived_effective_D4', dd['effective_d4_kcal_mol_from_full_minus_electronic'], 'kcal/mol', 'archived_context'),
        ('transferred_AMOEBA', str(dd['original_amoeba_transfer_value_kcal_mol']), 'kcal/mol', 'archived_context'),
        ('model_minus_archived_full', dd['original_amoeba_minus_published_reference_kcal_mol'], 'kcal/mol', 'archived_context'),
    ]
    for name in ['new_full_quantum', 'CT_suppressed_quantum', 'CT', 'polarization', 'combined_frozen',
                 'new_fragment_reference_residual', 'new_model_minus_full_quantum']:
        table.append((name, '', 'kcal/mol', 'unavailable_pilot_incomplete'))
    for name in ['native_energy_reconstruction_residual', 'new_full_vs_archive_residual']:
        table.append((name, '', 'hartree', 'unavailable_pilot_incomplete'))
    table.append(('CT_numerical_uncertainty', '', 'kcal/mol', 'not_established'))
    with (D / 'RESULTS_TABLE.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['quantity', 'value', 'units', 'status'])
        writer.writerows(table)
    result = {
        'status': 'bounded_diagnostic_inconclusive_pilot_runtime_limit', 'created_utc': stamp,
        'scientific_question_answered': False, 'pilot_admitted': False, 'remaining_three_released': False,
        'new_admitted_scientific_cells': 0, 'geometry_inventory': scope,
        'cause': 'First charged-La fragment did not converge before the declared native process cap; no ALMO chain completed.',
        'not_established': ['CT magnitude or sign', 'negligible CT', 'force-field missing mechanism',
                            'fundamental CP2K or Hamiltonian incompatibility', 'solution or protein selectivity'],
        'native_terminal_markers': markers, 'native_first_fragment_primary_steps': len(ot),
        'last_primary_step': ot[-1], 'printed_last_residual': '0.00000040',
        'fragment_EPS_SCF': '1e-10', 'native_exit_code': native['returncode'],
        'native_timed_out': True, 'native_cap_seconds': native['cap_seconds'],
        'native_elapsed_seconds_including_termination': native['elapsed_seconds'],
        'parser_result_absent_reason': 'The timeout exception precedes parser invocation; no missing success result is inferred.',
        'run_result': {'path': str((dest / 'RUN_RESULT.json').relative_to(D)), 'sha256': sha(dest / 'RUN_RESULT.json')},
        'mirrored_raw_artifacts': mirrored, 'config_sha256': sha(CONFIG), 'verified_execution_pins': pins,
        'goal_sha256': sha(GOAL), 'energy_ledger_sha256': sha(LEDGER),
        'table_sha256': sha(D / 'RESULTS_TABLE.csv'), 'costs_sha256': sha(D / 'COSTS.json'),
        'recommendation': 'Stop this attempt as inconclusive; do not use it to justify force-field changes or automatically extend the CP2K campaign.',
        'further_experiments_authorized_by_this_result': False,
    }
    save(D / 'FINAL_RESULT.json', result)
    save(TASK / 'ROOT_REVIEW.json', {'reviewed_utc': stamp, 'job_id': '21623', 'scheduler_state': 'FAILED',
         'cause': result['cause'], 'actual_callback_delivery_claimed': False,
         'final_result_path': str(D / 'FINAL_RESULT.json'), 'final_result_sha256': sha(D / 'FINAL_RESULT.json'),
         'action': 'Bounded test ended inconclusive; retain all attempts and leave remaining three unrun.'})
    print(json.dumps({'final_status': result['status'], 'allocated_cpu_seconds': allocated,
                      'new_admitted_cells': 0, 'verified_pins': len(pins)}, indent=2))


if __name__ == '__main__':
    main()
