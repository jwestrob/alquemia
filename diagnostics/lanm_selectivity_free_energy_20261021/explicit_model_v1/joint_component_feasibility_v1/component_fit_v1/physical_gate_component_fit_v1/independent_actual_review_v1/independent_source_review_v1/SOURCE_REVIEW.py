#!/usr/bin/env python3
"""Finite source/metadata receipt only; never import the reviewed reviewer."""
import ast
import datetime
import hashlib
import json
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
checks = []
pins = {}

def check(ok, label):
    if not ok:
        raise ValueError(label)
    checks.append(label)

def pin(path):
    path = Path(path).resolve()
    value = {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    pins[str(path)] = value
    return value

def read(path):
    return json.loads(Path(path).read_text())

def verify_items(value):
    # Follow only literal pin objects in the given source manifests, not files they reference.
    if isinstance(value, dict):
        if set(value) == {'path', 'sha256'}:
            check(pin(value['path']) == value, 'unchanged pin: ' + value['path'])
        else:
            for child in value.values():
                verify_items(child)
    elif isinstance(value, list):
        for child in value:
            verify_items(child)

def functions(path):
    return {n.name: ast.dump(n, include_attributes=False)
            for n in ast.parse(path.read_text()).body if isinstance(n, ast.FunctionDef)}

def main():
    wall, cpu = time.monotonic(), time.process_time()
    ready = read(PACKAGE / 'READY_v2.json')
    old_manifest = read(PACKAGE / 'MANIFEST.json')
    new_manifest = read(PACKAGE / 'MANIFEST_v2.json')
    for obj in (ready, old_manifest, new_manifest):
        verify_items(obj)
    check(ready['reviewer_implementation']['sha256'] == 'c1cb825804b19c3700e280ea328568a9d9f35882adee1ab26901a270b01fb268', 'frozen v2 source')
    check(ready['terminal_binding'] is None and ready['actual_review_executed'] is False, 'terminal remains unbound')
    check(ready['argv'][1] == ready['reviewer_implementation']['path'], 'argv uses v2')
    i = ready['argv'].index('--derived-replay-helper')
    check(ready['argv'][i + 1] == ready['derived_replay_helper']['path'], 'exact qualified arithmetic-helper argv')
    check(all(v == '1' for k, v in ready['environment'].items() if k != 'PYTHONDONTWRITEBYTECODE'), 'one-thread declared environment')
    old = functions(PACKAGE / 'REVIEW.py')
    new = functions(PACKAGE / 'REVIEW_v2.py')
    check(set(new) - set(old) == {'check_proposal_counts'}, 'one new counter helper')
    check([name for name in old if old[name] != new[name]] == ['run'], 'only run changes among inherited functions')
    inherited = read(ready['controls']['path'])
    changed = read(ready['changed_policy_controls']['path'])
    check(inherited['status'] == 'passed_artificial_terminal_review_controls' and inherited['error'] is None and len(inherited['checks']) == 15, '15 existing controls read only')
    check(changed['status'] == 'passed_proposal_count_source_controls' and changed['error'] is None and len(changed['checks']) == 10, '10 existing counter controls read only')
    check(changed['source'] == ready['reviewer_implementation'], 'counter controls bind final v2')
    regression = PACKAGE.parent / 'postselection_regression_v1' / 'RUN.py'
    text = regression.read_text()
    for literal in ('passed_physical_gate_component_fit_actual_review', 'selected_native_admission_passed', 'configuration', 'implementation', 'selected'):
        check(literal in text, 'postselection requires ' + literal)
    report = pin(HERE / 'REPORT.md')
    result = {
        'status': 'passed_physical_gate_terminal_reviewer_source_review',
        'implementation': ready['reviewer_implementation'],
        'reviewed_ready': pin(PACKAGE / 'READY_v2.json'),
        'previous_ready': ready['previous_ready'],
        'previous_implementation': old_manifest['sources']['REVIEW.py'],
        'manifest': ready['manifest'],
        'contract': ready['contract'],
        'derived_replay_helper': ready['derived_replay_helper'],
        'campaign_ready': ready['campaign_ready'],
        'inherited_controls': ready['controls'],
        'changed_policy_controls': ready['changed_policy_controls'],
        'postselection_implementation': pin(regression),
        'source_review_implementation': pin(__file__),
        'report': report,
        'reviewed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'findings': {
            'journal_and_new_receipt_accounting': 'passed source review',
            'complete_batch_selection_and_tie_policy': 'passed source review',
            'selected_native_admission_and_original25': 'passed source review',
            'terminal_failure_nonadmission': 'passed source review',
            'logical_versus_real_proposal_counts': 'resolved by immutable v2',
            'postselection_status_contract': 'compatible; actual binding remains required'
        },
        'required_actual_pass_status': ready['expected_status'],
        'required_actual_fields': ready['required_contract'],
        'checks': checks,
        'unique_verified_pin_count': len(pins),
        'source_only': True,
        'campaign_progress_or_outcomes_read': False,
        'controls_rerun': False,
        'reviewer_executed': False,
        'new_native_calls': 0,
        'new_optimizer_calls': 0,
        'source_authorship_disclosure': 'Planck reviewed the execution source and authored the downstream trimer adapter; Singer authored this terminal reviewer and driver/adapter. CORE has its separate review.',
        'limits': [
            'Source readiness only; no terminal campaign or selected candidate is admitted by this receipt.',
            'Actual execution must retain the exact READY_v2 argv and its pinned derived arithmetic helper.',
            'Historical native evidence is reused; all new receipts are audited but raw scientific replay is limited to newly generated final-selected cells.',
            'Failure audit retains cost and invocation evidence without granting selected-native admission.',
            'Measured login CPU/wall values are reported, not independently reconstructed scheduler allocation or physical validation.'
        ],
        'cost': {'wall_seconds': time.monotonic() - wall, 'CPU_seconds': time.process_time() - cpu, 'threads': 1}
    }
    out = HERE / 'REVIEW.json'
    with out.open('x') as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'review': pin(out), 'checks': len(checks), 'unique_pins': result['unique_verified_pin_count'], 'cost': result['cost']}))

if __name__ == '__main__':
    main()
