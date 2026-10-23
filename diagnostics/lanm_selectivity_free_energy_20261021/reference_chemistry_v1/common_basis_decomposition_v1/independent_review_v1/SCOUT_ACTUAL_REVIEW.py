#!/usr/bin/env python3
"""Saved-only independent admission of the two declared common-basis scouts."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import time

HERE = Path(__file__).resolve().parent
BRANCH = HERE.parent
ROOT = BRANCH.parents[3]
WORK = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/common_basis_decomposition_v1'


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    tick, cpu = time.monotonic(), time.process_time()
    review = json.loads((HERE / 'REVIEW.json').read_text())
    assert review['implementation'] == pin(BRANCH / 'REFERENCE.py')
    spec = importlib.util.spec_from_file_location('scout_reference_admission', BRANCH / 'REFERENCE.py')
    a = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(a)
    panel_path = WORK / 'prepared_v1/PANEL.json'
    panel = a.load(panel_path)
    checks = {}

    def check(label, condition):
        assert condition, label
        checks[label] = True

    paths = {cid: WORK / 'run_v1/tasks' / cid / 'execution/RESULT.json' for cid in a.SCOUTS}
    # Read neither live trajectory nor a partial native output.
    check('both_terminal_result_files_exist', all(path.is_file() for path in paths.values()))
    raw = {cid: json.loads(path.read_text()) for cid, path in paths.items()}
    check('both_terminal_statuses_complete', all(r['status'] == 'complete' for r in raw.values()))
    rows = {}
    for cid, path in paths.items():
        r = raw[cid]
        cell = next(c for c in panel['cells'] if c['cell_id'] == cid)
        admitted = a.admit_result(path, cell, panel, panel_path, require_scout=True)
        text = a.b.verify(r['output']).read_text()
        start = a.read(r['start'])
        role = cell['role']
        electrons, charge = (18, 3) if role == 'metal_water' else (10, 0)
        prefix = role + ':'
        check(prefix + 'all_declared_native_checks', len(r['checks']) == 19 and all(r['checks'].values()))
        for name, pattern, expected in [
            ('electrons', r'Number of Electrons\s+NEL\s+\.+\s+(\d+)', electrons),
            ('charge', r'Total Charge\s+Charge\s+\.+\s+(-?\d+)', charge),
            ('multiplicity', r'Multiplicity\s+Mult\s+\.+\s+(\d+)', 1),
            ('orbital_functions', r'Number of basis functions\s+\.+\s+(\d+)', 434),
            ('aux_header', r'# of basis functions in Aux-J\s+\.+\s+(\d+)', 1008),
            ('aux_integral', r'Number of AuxJ basis functions\s+\.+\s+(\d+)', 1008),
        ]:
            values = [int(v) for v in re.findall(pattern, text)]
            check(prefix + name, bool(values) and set(values) == {expected})
        finals = re.findall(r'^FINAL SINGLE POINT ENERGY\s+([-+0-9.]+)\s*$', text, re.M)
        totals = re.findall(r'^Total Energy\s+:\s+([-+0-9.]+)\s+Eh', text, re.M)
        check(prefix + 'raw_final_energy', len(finals) == 1 and float(finals[0]) == r['energy_hartree'])
        check(prefix + 'raw_electronic_energy', len(totals) == 1 and float(totals[0]) == r['electronic_energy_hartree'])
        check(prefix + 'electronic_closure', abs(r['energy_hartree'] - r['electronic_energy_hartree']) < 1e-10 and r['D4_energy_hartree'] == 0.)
        check(prefix + 'real_ghost_nuclear_charges', [v['effective_Z'] for v in r['atomic_table']] == cell['expected_effective_nuclear_charges'])
        check(prefix + 'basis_each_center', r['basis_functions_by_center'] == [47] + [31, 6, 6] * 9)
        check(prefix + 'ecp_core', r['ecp_core_counts'] == ([55] if role == 'metal_water' else []))
        check(prefix + 'ecp_real_center_only', r['ecp_assigned_centers'] == ([0] if role == 'metal_water' else []))
        check(prefix + 'literal_28_center_coordinates', len(r['atomic_table']) == 28 and r['coordinate_max_error_A'] < 1e-6)
        check(prefix + 'effective_vs_physical_spin', cell['physical_multiplicity_of_parent'] == 6 and cell['multiplicity'] == 1)
        check(prefix + 'exact_input_rerender', a.b.verify(r['input']).read_text() == a.deck(panel, cell, 4, start['maxcore_MiB']))
        check(prefix + 'one_four_rank_call', r['native_program_invocations_attempted'] == 1 and start['ranks'] == 4)
        check(prefix + 'native_elapsed_finite', math.isfinite(r['elapsed_seconds']) and r['elapsed_seconds'] > 0 and r['allocated_cpu_seconds'] == 4 * r['elapsed_seconds'])
        criteria = {v['name']: v for v in r['printed_scf_criteria']}
        check(prefix + 'scf_metadata_retained', len(criteria) == 6 and bool(r['native_scf_stop_signals']))
        check(prefix + 'no_state_accuracy_claim', r['physical_state_qualified'] is False and r['method_accuracy_qualified'] is False)
        rows[cid] = {
            'result': pin(path), 'output': r['output'], 'input': r['input'], 'start': r['start'], 'release': r['release'],
            'job_id': r['job_id'], 'role': role, 'electrons': electrons, 'charge': charge,
            'orbital_functions': 434, 'auxiliary_functions': 1008,
            'electronic_energy_hartree': r['electronic_energy_hartree'], 'D4_energy_hartree': r['D4_energy_hartree'],
            'printed_scf_criteria': r['printed_scf_criteria'], 'native_scf_stop_signals': r['native_scf_stop_signals'],
            'final_COSX_correction_hartree': r['final_COSX_correction_hartree'],
            'total_energy_after_final_integration_hartree': r['total_energy_after_final_integration_hartree'],
            'native_warning_lines': [line.strip() for line in text.splitlines() if 'WARNING:' in line],
            'native_elapsed_seconds': r['elapsed_seconds'], 'native_elapsed_times_four_ranks_CPU_seconds': r['allocated_cpu_seconds'],
            'raw_admission_status': admitted['status'],
        }
    result = {
        'status': 'passed_independent_actual_scout_review', 'implementation': pin(__file__),
        'adapter': panel['adapter'], 'panel': pin(panel_path), 'independent_source_review': pin(HERE / 'REVIEW.json'),
        'scout_results': {cid: pin(path) for cid, path in paths.items()}, 'rows': rows,
        'checks': checks, 'passed_checks': len(checks), 'declared_native_scope': 68, 'completed_scouts': 2,
        'remaining_cells': 66, 'ready_for_root_remaining_release': True,
        'limitations': [
            'Both full-parent ghost masks are supported by these actual W9 Dy scouts; other cells still require their own unchanged native admission.',
            'Native stopping signals and all printed residual flags are retained. Native completion is not all-criterion, physical-state or basis-limit accuracy qualification.',
            'No energy decomposition, target change, force qualification or model promotion follows from two scouts.',
            'Native elapsed times ranks is an allocation-time proxy, not measured CPU utilization; root retains scheduler accounting separately.'
        ],
        'cost': {'review_wall_seconds': time.monotonic() - tick, 'review_parent_CPU_seconds': time.process_time() - cpu,
                 'new_native_calls': 0, 'local_threads': 1,
                 'reviewed_native_elapsed_times_ranks_CPU_seconds': sum(v['native_elapsed_times_four_ranks_CPU_seconds'] for v in rows.values())},
    }
    target = HERE / 'SCOUT_ACTUAL_REVIEW.json'
    assert not target.exists()
    target.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': result['status'], 'review': pin(target), 'checks': len(checks), 'new_native_calls': 0}))


if __name__ == '__main__':
    main()
