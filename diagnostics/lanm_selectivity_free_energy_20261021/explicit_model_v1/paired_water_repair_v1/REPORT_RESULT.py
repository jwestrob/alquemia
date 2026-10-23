#!/usr/bin/env python3
"""Review/render a terminal pilot result using saved fields and hashes only."""
import argparse
import hashlib
import json
from pathlib import Path
import time


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def main(result_path, output, completion=None):
    started = time.monotonic()
    p = Path(result_path).resolve()
    result = json.loads(p.read_text())
    if result['status'] not in ('finite_pilot_complete', 'failed_or_unavailable'):
        print(json.dumps({'status': 'not_terminal', 'producer_status': result['status'],
                          'search_status': result.get('search', {}).get('status'), 'no_report_written': True}))
        return 2
    verified = {}
    def check(record):
        key = record['path']
        if key not in verified:
            verified[key] = pin(key)
        if verified[key] != record:
            raise ValueError('Changed evidence: '+key)
        return Path(key)
    config = json.loads(check(result['config']).read_text())
    for key in ('implementation', 'specification', 'objective', 'search', 'evaluator',
                'evaluator_config', 'body_helper', 'water48', 'root_release', 'preparation_result'):
        check(config[key])
    if result['implementation'] != config['implementation']:
        raise ValueError('Driver source association differs')
    evaluator = result.get('evaluator', {})
    ledger = json.loads(check(evaluator['ledger']).read_text()) if evaluator else None
    if ledger and ledger['native_program_invocations'] != result['new_native_program_invocations']:
        raise ValueError('Native counter differs')
    if result['new_native_program_invocations'] > 1920:
        raise ValueError('Native ceiling exceeded')
    search = result.get('search', {})
    evaluations = search.get('evaluations', {})
    native_records = {}
    for evaluation in evaluations.values():
        for gid, row in (evaluation.get('model') or {}).items():
            provenance = row['provenance']
            authority = provenance['native_authority']
            check(authority)
            for name in ('analyze_output', 'gradient_output'):
                check(row['native'][name])
            for v in provenance['source_input'].values():
                check(v)
            if provenance['reused_baseline']:
                continue
            if authority['path'] not in native_records:
                native = json.loads(Path(authority['path']).read_text())
                if native['status'] != 'complete' or native['cell_id'] != gid:
                    raise ValueError('Native authority unavailable or wrong cell')
                if native['energy'] != row['energy'] or native['forces'] != row['forces']:
                    raise ValueError('Saved native arrays differ')
                for v in native['input'].values():
                    check(v)
                if len(native['receipts']) != 2 or any(v['returncode'] != 0 for v in native['receipts']):
                    raise ValueError('Native invocation receipt differs')
                native_records[authority['path']] = native
    complete = result['status'] == 'finite_pilot_complete'
    if complete:
        if len(evaluations) != search['distinct_vectors'] or len(evaluations) > 241:
            raise ValueError('Vector denominator differs')
        if any(row['status'] != 'complete' for row in evaluations.values()):
            raise ValueError('Complete pilot contains unavailable vectors')
        if any(not poll['complete'] for poll in search['polls']):
            raise ValueError('Complete pilot contains a partial selected poll')
        if 2*len(native_records) != result['new_native_program_invocations']:
            raise ValueError('Unique native receipts do not close the call counter')
        best = evaluations[search['best_key']]
        if best['parameters'] != result['candidate_parameters'] or best['score'] != result['selected_score']:
            raise ValueError('Selected candidate differs')
        if len(result['selected_model']) != 4 or result['selected_model'] != best['model']:
            raise ValueError('Selected four-cell denominator differs')
        for entry in result['candidate_keys'].values():
            check(entry['key'])
    scheduler = None
    if completion:
        scheduler = json.loads(Path(completion).read_text())
    summary = {'status': 'saved_result_review_passed', 'implementation': pin(__file__), 'result': pin(p),
        'producer_status': result['status'], 'producer_error': result.get('error'),
        'verified_unique_files': len(verified), 'new_native_program_invocations': result['new_native_program_invocations'],
        'unique_complete_native_cells': len(native_records), 'distinct_vectors': len(evaluations),
        'complete_polls': len(search.get('polls', [])), 'search_stopping_rule': search.get('search_stopping_rule'),
        'numerical_optimization_complete': result.get('numerical_optimization_complete'),
        'optimizer_budget_limited': result.get('optimizer_budget_limited'),
        'development_checks_passed': result.get('development_checks_passed'),
        'identifiability_and_bounds_passed': result.get('identifiability_and_bounds_passed'),
        'near_bound_by_parameter': search.get('near_bound_by_parameter'), 'sensitivity': search.get('sensitivity'),
        'baseline_score': result.get('baseline_score'), 'selected_score': result.get('selected_score'),
        'selected_method_sensitivity': result.get('selected_method_sensitivity'),
        'cost': {k: result.get(k) for k in ('elapsed_seconds', 'parent_process_CPU_seconds',
            'aggregate_native_child_user_CPU_seconds', 'aggregate_native_child_system_CPU_seconds')},
        'scheduler_snapshot': {'pin': pin(completion), 'status': scheduler.get('status'),
            'records': scheduler.get('records')} if scheduler else None,
        'new_scientific_calculations': 0, 'physical_accuracy_qualified': False,
        'validation_outcomes_opened': False, 'automatic_extension_or_validation_release': False}
    lines = ['# Paired water pilot: saved result review', '',
        f"Producer status: **{result['status']}**. Saved source/native hash and receipt review passed.", '',
        f"- Distinct parameter vectors: {len(evaluations)}/241; complete polls: {len(search.get('polls', []))}.",
        f"- Search stop: `{search.get('search_stopping_rule')}`; local stopping rule met: {result.get('numerical_optimization_complete')}.",
        f"- New native invocations: {result['new_native_program_invocations']}/1920, covering {len(native_records)} unique complete cells.",
        f"- Recorded wall time: {result.get('elapsed_seconds')} s. Parent CPU: {result.get('parent_process_CPU_seconds')} s.",
        f"- Aggregate native child CPU: user {result.get('aggregate_native_child_user_CPU_seconds')} s; system {result.get('aggregate_native_child_system_CPU_seconds')} s. Overlapping worker receipts are not summed.", '']
    if complete:
        old, selected = result['baseline_score'], result['selected_score']
        lines += ['## Separate development controls', '',
            '| Response | Initial RMS | Selected RMS | Selected maximum component | Pass |',
            '|---|---:|---:|---:|---|']
        for kind in ('La', 'Dy', 'paired'):
            a, b = old['controls'][kind], selected['controls'][kind]
            lines.append(f"| {kind} | {a['category_balanced_rigid_RMS_kcal_mol_A']:.6f} | {b['category_balanced_rigid_RMS_kcal_mol_A']:.6f} | {b['max_body_component_error_kcal_mol_A']:.6f} | {b['passed']} |")
        lines += ['', 'Force units are kcal/mol/Å; torques are divided by1 Å.', '',
            f"Balanced exchange error: {old['balanced_exchange_error_kcal_mol']:.6f} → {selected['balanced_exchange_error_kcal_mol']:.6f} kcal/mol (threshold1).",
            f"Objective: {old['objective']:.6f} → {selected['objective']:.6f}. All separate controls pass: **{result['development_checks_passed']}**.", '',
            '## Parameters and method sensitivity', '',
            f"Selected parameters: `{result['candidate_parameters']}` (La radius/damping, Dy radius/damping).",
            f"Near-bound flags: `{search['near_bound_by_parameter']}`. Jacobian rank: {search['sensitivity']['rank']}/4; condition: {search['sensitivity']['condition']}. Weakly identified: {search['sensitivity']['weakly_identified']}.", '',
            '| Response | Selected dE/dq | UnCP PBE0-D4 residual | UnCP CC residual | CP CC residual |',
            '|---|---:|---:|---:|---:|']
        for kind, row in result['selected_method_sensitivity']['results'].items():
            comparisons = row['comparisons']
            lines.append(f"| {kind} | {row['model_dE_dq_kcal_mol_A']:.6f} | {comparisons['PBE0_D4_uncorrected']['residual_kcal_mol_A']:.6f} | {comparisons['CC_uncorrected']['residual_kcal_mol_A']:.6f} | {comparisons['CC_CP']['residual_kcal_mol_A']:.6f} |")
        lines += ['', 'Positive q moves rigid water outward; derivatives and residuals have kcal/mol/Å units.',
            f"CC residual deterioration flag: **{result['selected_method_sensitivity']['CC_residual_deterioration_flag']}**. Full sign and method comparisons remain in REVIEW.json."]
    else:
        lines += [f"Producer error: `{result.get('error')}`. Partial records remain unavailable; no replacement candidate is selected."]
    lines += ['', '## Interpretation limits', '',
        'A complete run is distinct from passing the development controls. A budget stop does not establish converged optimization or failure of the parameter family. A local step stop does not establish a global optimum. Rank/bound flags concern this small consumed design, not statistical parameter uncertainty.', '',
        'This is representability against one uncorrected PBE0 surface. CC/CP comparisons remain separate method sensitivity; neither bulk hydration, protein donor transfer, explicit-f state accuracy nor affinity is qualified. Reserved validation stays closed. No retry, extension, parameter change or production trajectory is released.', '',
        f"Result: `{p}`", f"SHA256: `{summary['result']['sha256']}`", '']
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    summary['review_elapsed_seconds'] = time.monotonic()-started
    (output/'REVIEW.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    (output/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({'review': pin(output/'REVIEW.json'), 'report': pin(output/'REPORT.md')}))
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--result', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--completion-event')
    args = parser.parse_args()
    raise SystemExit(main(args.result, args.output, args.completion_event))
