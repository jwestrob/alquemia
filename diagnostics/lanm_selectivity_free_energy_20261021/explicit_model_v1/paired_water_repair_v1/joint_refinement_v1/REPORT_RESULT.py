#!/usr/bin/env python3
"""Saved-only joint-refinement report: hashes, receipts, raw parsing and scalars."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import time

CELLS = {'water_central_La', 'water_central_Dy', 'aquo8_central_La', 'aquo8_central_Dy'}


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def token(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def brief(row):
    score = row['score']
    return {k: score[k] for k in ('objective', 'balanced_exchange_error_kcal_mol', 'controls', 'development_checks_passed')}


def main(path, output, completion=None):
    started = time.monotonic()
    source = Path(path).resolve(); result = json.loads(source.read_text())
    if result['status'] != 'bounded_joint_refinement_complete':
        print(json.dumps({'status': 'not_completed', 'producer_status': result['status'], 'no_report_written': True})); return 2
    verified = {}
    def check(p):
        if p['path'] not in verified: verified[p['path']] = pin(p['path'])
        if verified[p['path']] != p: raise ValueError('Changed evidence: ' + p['path'])
        return Path(p['path'])
    def read(p): return json.loads(check(p).read_text())
    cfg = read(result['config'])
    for p in cfg.values():
        if isinstance(p, dict) and set(p) == {'path', 'sha256'}: check(p)
    assert result['implementation'] == cfg['implementation']
    ecfg = read(cfg['evaluator_config'])
    helper_pin = ecfg['native_helper']; check(helper_pin)
    spec = importlib.util.spec_from_file_location('saved_native_parser', helper_pin['path'])
    helper = importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
    ledger = read(result['evaluator']['ledger'])
    assert ledger['native_program_invocations'] == result['new_native_program_invocations'] <= 1152
    for parent in result['protected_parents']: check(parent)
    seed_paths = set()
    for seed in result['seeded_cache']:
        check(seed['authority']); seed_paths.add(seed['authority']['path'])
        assert Path(seed['symlink']).is_symlink()
        assert Path(seed['symlink']).resolve() == Path(seed['authority']['path']).parent
    vectors = result['vectors']; native_records = {}; raw_parsed = set()
    for key, row in vectors.items():
        assert row['status'] == 'complete' and token(row['parameters']) == key and set(row['model']) == CELLS
        for gid, cell in row['model'].items():
            authority = cell['provenance']['native_authority']; check(authority)
            for p in cell['provenance']['source_input'].values(): check(p)
            for name in ('analyze_output', 'gradient_output'): check(cell['native'][name])
            assert cell['energy'] == cell['native']['total_energy_kcal_mol']
            assert cell['forces'] == cell['native']['force_kcal_mol_A']
            if not cell['provenance']['reused_baseline']:
                raw = read(authority)
                assert raw['status'] == 'complete' and raw['cell_id'] == gid
                assert raw['energy'] == cell['energy'] and raw['forces'] == cell['forces'] and raw['native'] == cell['native']
                for p in raw['input'].values(): check(p)
                assert len(raw['receipts']) == 2 and all(r['returncode'] == 0 for r in raw['receipts'])
                for receipt in raw['receipts']: check(receipt['output'])
                native_records[authority['path']] = raw
    new_native = {p: r for p, r in native_records.items() if p not in seed_paths}
    assert len(new_native) * 2 == result['new_native_program_invocations']
    for p in new_native:
        assert Path(p).is_relative_to(source.parent/'cache/cells')
    assert len(vectors)-1 == result['new_parameter_vectors'] <= 144
    assert result['completed_outer_steps'] == len(result['iterations']) <= 12
    assert result['original_results_and_ledgers_unchanged'] and result['candidate_frozen']
    selected = result['selected']; assert selected == vectors[result['current_key']]
    selected_native = {}
    for gid, row in selected['model'].items():
        raw = helper.parse(check(row['native']['analyze_output']), check(row['native']['gradient_output']),
                           {'atoms': len(row['forces']), 'derivative_test': False})
        assert raw == row['native']
        raw_parsed.add(gid)
        authority = row['provenance']['native_authority']; original = read(authority)
        selected_native[gid] = {'authority': authority, 'input': original['input'],
            'analyze_output': raw['analyze_output'], 'gradient_output': raw['gradient_output'],
            'energy_kcal_mol': row['energy'], 'atoms': len(row['forces']), 'raw_replay_exact': True}
    iterations = []; all_passing = []
    current = token(result['initial']['parameters'])
    for entry in result['iterations']:
        assert entry['status'] == 'complete' and entry['origin_key'] == current
        assert len(entry['derivative_keys']) == 8 and len(entry['trial_keys']) == 4
        assert len(entry['derivative_slots']) == 8 and len(entry['trial_slots']) == 4
        assert all(k in vectors for k in entry['derivative_keys'] + entry['trial_keys'])
        best = current
        for k in entry['trial_keys']:
            if vectors[k]['score']['objective'] < vectors[best]['score']['objective']: best = k
        assert best == entry['selected_key']
        trials = []
        for slot, k in zip(entry['trial_slots'], entry['trial_keys']):
            row = vectors[k]; assert slot['parameters'] == row['parameters']
            rec = dict(id=slot['id'], key=k, selected=(k == best), **brief(row))
            trials.append(rec)
            if row['score']['development_checks_passed'] and k != best:
                all_passing.append({'outer_step': entry['outer_step'], **rec})
        singular = entry['linear_directions']['singular_values']
        iterations.append({'outer_step': entry['outer_step'], 'origin_key': current,
            'origin': brief(vectors[current]), 'selected_key': best, 'selected': brief(vectors[best]),
            'derivative_widths': entry['derivative_actual_widths'], 'trials': trials,
            'origin_Jacobian_rank': entry['linear_directions']['rank'], 'origin_singular_values': singular,
            'origin_condition': singular[0]/singular[-1] if singular[-1] else None})
        current = best
    assert current == result['current_key']
    assert result['development_checks_passed'] == selected['score']['development_checks_passed']
    if result['stopping_reason'] == 'all_existing_controls_passed': assert result['development_checks_passed']
    params = []
    for j, p in enumerate(result['family']['parameters']):
        x = selected['parameters'][j]; z = selected['normalized_parameters'][j]
        assert p['lower'] <= x <= p['upper']
        params.append({'name': p['name'], 'initial_scout': result['initial']['parameters'][j], 'selected': x,
            'move_from_scout': x-result['initial']['parameters'][j], 'original_family_initial': p['initial'],
            'lower': p['lower'], 'upper': p['upper'], 'normalized_selected': z,
            'normalized_distance_to_nearest_bound': min(z, 1-z), 'exact_bound_hit': x in (p['lower'], p['upper'])})
    event = json.loads(Path(completion).read_text()) if completion else None
    summary = {'status': 'saved_result_review_passed', 'reporter': pin(__file__), 'result': pin(source),
        'config': result['config'], 'implementation': result['implementation'], 'producer_status': result['status'],
        'stopping_reason': result['stopping_reason'], 'initial': brief(result['initial']), 'selected': brief(selected),
        'selected_parameters': selected['parameters'], 'parameters': params, 'fixed_pair_epsilon': result['family']['fixed_pair_epsilon'],
        'iterations': iterations, 'passing_unselected_trials': all_passing,
        'all_selected_native_raw_replays_exact': len(raw_parsed) == 4, 'selected_native': selected_native,
        'verified_unique_files': len(verified), 'inherited_native_cells_seeded': len(seed_paths),
        'new_native_cells': len(new_native), 'new_native_program_invocations': result['new_native_program_invocations'],
        'new_parameter_vectors': result['new_parameter_vectors'], 'complete_outer_steps': len(iterations),
        'method_sensitivity': result['selected_method_sensitivity'],
        'cost': {k: result[k] for k in ('elapsed_seconds', 'parent_process_CPU_seconds',
            'aggregate_native_child_user_CPU_seconds', 'aggregate_native_child_system_CPU_seconds')},
        'scheduler_records': event.get('records') if event else None,
        'final_point_Jacobian_not_evaluated': True, 'global_optimum_claimed': False,
        'new_model_evaluations': 0, 'physical_accuracy_qualified': False, 'validation_outcomes_opened': False}
    lines = ['# Bounded joint water refinement: actual result', '',
        f"All separate development controls pass. Stop: `{result['stopping_reason']}` after {len(iterations)} complete outer steps.", '',
        f"Objective {summary['initial']['objective']:.9f} → {summary['selected']['objective']:.9f}; balanced exchange error {summary['selected']['balanced_exchange_error_kcal_mol']:.9f} kcal/mol (limit 1).", '',
        '| Control | Selected RMS | Selected maximum | Pass |', '|---|---:|---:|---|']
    for name, row in summary['selected']['controls'].items():
        lines.append(f"| {name} | {row['category_balanced_rigid_RMS_kcal_mol_A']:.6f} | {row['max_body_component_error_kcal_mol_A']:.6f} | {row['passed']} |")
    lines += ['', 'Rigid-response units are kcal/mol/Å; torques use the unchanged 1 Å divisor. RMS/max limits remain 2/5.', '',
        '## Completed steps', '', '| Step | Selected trial | Objective | La max | Dy max | Paired max | All pass |', '|---|---|---:|---:|---:|---:|---|']
    for entry in iterations:
        sc = entry['selected']; names = [t['id'] for t in entry['trials'] if t['selected']]
        lines.append(f"| {entry['outer_step']} | {','.join(names) or 'current'} | {sc['objective']:.9f} | {sc['controls']['La']['max_body_component_error_kcal_mol_A']:.6f} | {sc['controls']['Dy']['max_body_component_error_kcal_mol_A']:.6f} | {sc['controls']['paired']['max_body_component_error_kcal_mol_A']:.6f} | {sc['development_checks_passed']} |")
    lines += ['', f"Passing but unselected trial slots: {len(all_passing)}. Every trial's controls and all such rows are retained in REVIEW.json. Selection follows the original objective with current-first ties.", '',
        '## Parameters and cost', '', '| Parameter | Scout start | Selected | Move | Bound hit |', '|---|---:|---:|---:|---|']
    for p in params: lines.append(f"| {p['name']} | {p['initial_scout']:.9f} | {p['selected']:.9f} | {p['move_from_scout']:+.9f} | {p['exact_bound_hit']} |")
    lines += ['', f"New vectors: {result['new_parameter_vectors']}/144; native invocations: {result['new_native_program_invocations']}/1152 across {len(new_native)} new cells; {len(seed_paths)} inherited cells preserved.",
        f"Recorded wall {result['elapsed_seconds']:.3f} s; parent CPU {result['parent_process_CPU_seconds']:.3f} s; aggregate native child user/system CPU {result['aggregate_native_child_user_CPU_seconds']:.3f}/{result['aggregate_native_child_system_CPU_seconds']:.3f} s. Concurrent receipt CPU deltas are not summed.", '',
        'All selected raw printed energy/gradient records replay exactly; sources, inputs, native receipts and parent ledgers retain their pins.', '',
        '## Method sensitivity and limits', '', '| Response | dE/dq | UnCP PBE0 residual | UnCP CC residual | CP CC residual |', '|---|---:|---:|---:|---:|']
    for name, row in result['selected_method_sensitivity']['results'].items():
        c = row['comparisons']; lines.append(f"| {name} | {row['model_dE_dq_kcal_mol_A']:.6f} | {c['PBE0_D4_uncorrected']['residual_kcal_mol_A']:.6f} | {c['CC_uncorrected']['residual_kcal_mol_A']:.6f} | {c['CC_CP']['residual_kcal_mol_A']:.6f} |")
    lines += ['', 'Positive q moves water outward; units are kcal/mol/Å. CC/CP comparisons remain separate from the objective. The CC deterioration flag is false relative to the original model, while individual CC residuals remain substantial.', '',
        'This is a controls-pass stop on four consumed water fixtures, not global optimality or physical qualification. Both stored Jacobians have rank 4 at their iteration origins; the final selected point has no newly evaluated Jacobian. Protein/donor transfer, bulk hydration, electronic-state accuracy and affinity remain unqualified. Reserved validation is unopened; this report releases no further evaluation.', '',
        f"Result SHA256: `{summary['result']['sha256']}`", '']
    out = Path(output); out.mkdir(parents=True, exist_ok=False)
    if event: (out/'COMPLETION_SNAPSHOT.json').write_text(json.dumps(event, indent=2)+'\n')
    summary['saved_review_elapsed_seconds'] = time.monotonic()-started
    (out/'REVIEW.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    (out/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({'review': pin(out/'REVIEW.json'), 'report': pin(out/'REPORT.md'), 'result': summary['result']}))
    return 0


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--result', required=True); p.add_argument('--output', required=True); p.add_argument('--completion-event')
    a = p.parse_args(); raise SystemExit(main(a.result, a.output, a.completion_event))
