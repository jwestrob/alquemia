#!/usr/bin/env python3
"""One saved-data summary, restricted to the five already observed accepted steps."""
import datetime
import hashlib
import json
import math
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / 'workspaces').is_dir())
REL = Path('explicit_model_v1/joint_component_feasibility_v1/component_fit_v1')
D = ROOT / 'diagnostics/lanm_selectivity_free_energy_20261021' / REL
W = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021' / REL

def pin(p, raw=None):
    return {'path': str(p.resolve()), 'sha256': hashlib.sha256(p.read_bytes() if raw is None else raw).hexdigest()}

def read(p):
    raw = p.read_bytes()
    return json.loads(raw), pin(p, raw)

def summarize(m):
    gates = m['original25']; scores = m['original_block_scores']
    e = {k: v['error_kcal_mol'] for k, v in gates.items() if k.startswith('energy/')}
    f = {k: v for k, v in gates.items() if k.startswith('force/')}
    return {'key': m['key'], 'q': m['q'], 'component_RMS_kcal_mol': m['component_RMS'],
        'component_errors_kcal_mol': m['component_errors_kcal_mol'],
        'derived_C_errors_kcal_mol': {str(n): m['component_errors_kcal_mol']['P'+str(n)] + m['component_errors_kcal_mol']['B'+str(n)] for n in (8, 9)},
        'original_block_scores': scores,
        'electrical_RMS': {a: math.sqrt(scores[b]) for a, b in [('pair_Debye','pair_mu'), ('C_Debye','C_mu'), ('chi_a0_cubed','chi')]},
        'energy_errors_kcal_mol': e, 'force_gates': f,
        'maximum_absolute_energy_error_kcal_mol': max(abs(x) for x in e.values()),
        'maximum_rigid_RMS_kcal_mol_A': max(x['rigid_RMS_kcal_mol_A'] for x in f.values()),
        'maximum_body_component_error_kcal_mol_A': max(x['maximum_body_component_error_kcal_mol_A'] for x in f.values()),
        'producer_reported_physical_gates_passed': sum(x['passed'] for x in gates.values()),
        'producer_reported_full_development_target': m['full_development_target'],
        'producer_reported_substantive_recovery_milestone': m['substantive_recovery_milestone']}

def main():
    tick, cpu = time.monotonic(), time.process_time()
    current, current_pin = read(W / 'physical_gate_component_fit_v1/run_v1/execution/RESULT.json')
    prior, prior_pin = read(W / 'run_v1/execution/RESULT.json')
    trimer, trimer_pin = read(W / 'postselection_regression_v1/independent_actual_review_v1/REVIEW.json')
    accepted = []
    for it in current['iterations']:
        if it['outer_step'] > 5:
            continue
        assert it['status'] == 'complete'
        m = it['selection']['selected_metrics']; p = current['vectors'][m['key']]['artifact']
        assert pin(Path(p['path'])) == p
        accepted.append({'outer_step': it['outer_step'], 'artifact': p, 'metrics': summarize(m)})
    assert len(accepted) == 5 and accepted[-1]['metrics']['key'] == '2fb6edff18981d9d1210cfac16c6fde788fee8b4fdf4088c982ab939b41b9175'
    fields = ('candidate_energy_error', 'previous_selected_energy_error', 'baseline_energy_error', 'candidate_mu_error_norm_Debye', 'previous_selected_mu_error_norm_Debye', 'baseline_mu_error_norm_Debye')
    value = {'status': 'preliminary_interpretation_snapshot_not_terminal_admission',
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'implementation': pin(Path(__file__)),
        'live_source_observation': {'pin_at_read': current_pin, 'observed_status': current['status'], 'mutable_path_warning': 'This hash records bytes read once; the live path may advance. Only the already observed accepted steps 1 through 5 are summarized.'},
        'origin24': summarize(current['initial_metrics']), 'accepted_steps_1_through_5': accepted,
        'prior_terminal_result': prior_pin, 'prior_terminal': summarize(prior['selected_metrics']),
        'original_prefit_component_RMS_kcal_mol': current['initial_metrics']['original_initial_component_RMS'],
        'prior_trimer_review': trimer_pin, 'prior_trimer_status': trimer['status'],
        'prior_trimer_cycles': {k: {f: v[f] for f in fields} for k, v in trimer['cycles'].items()},
        'prior_trimer_paired_Dy_minus_La': trimer['paired_Dy_minus_La'],
        'prior_trimer_comparison_counts': trimer['comparison_counts'],
        'source_authorities': {k: pin(p) for k, p in {
            'new_declaration': D / 'physical_gate_component_fit_v1/ROOT_DECLARATION.json',
            'component_definitions': D.parent / 'PROTOCOL.md',
            'matched_conic_interpretation': D / 'ceiling_tradeoff_v1/interpretation_v1/REPORT.md',
            'matched_native_pilot_interpretation': D / 'ceiling_tradeoff_v1/nonlinear_check_v1/REPORT.md',
            'prior_trimer_report': W / 'postselection_regression_v1/independent_actual_review_v1/REPORT.txt',
            'terminal_reviewer_readiness': D / 'physical_gate_component_fit_v1/independent_actual_review_v1/READY_v2.json',
            'planned_trimer_readiness': D / 'physical_gate_component_fit_v1/postselection_regression_v1/READY.json'}.items()},
        'pending_bindings': {'terminal_result': None, 'terminal_review': None, 'postselection_result': None, 'postselection_review': None},
        'new_native_calls': 0, 'new_reference_calls': 0, 'new_optimizer_calls': 0, 'reserved_outcomes_read': False,
        'raw_native_readmission_performed': False,
        'cost': {'wall_seconds': time.monotonic() - tick, 'CPU_seconds': time.process_time() - cpu, 'threads': 1}}
    out = HERE / 'VALUES.json'
    with out.open('x') as handle:
        handle.write(json.dumps(value, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'values': pin(out), 'cost': value['cost'], 'observed_status': current['status']}))

if __name__ == '__main__':
    main()
