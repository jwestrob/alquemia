#!/usr/bin/env python3
"""Saved initial-ten arithmetic control; never reads the live fit result."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[6]
REL = Path('lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_coordination_joint_v1')
D = ROOT / 'diagnostics' / REL
W = ROOT / 'workspaces' / REL

def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def load_module(p, name):
    spec = importlib.util.spec_from_file_location(name, p)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def main():
    start, cpu = time.monotonic(), time.process_time()
    assert all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'))
    prep_path = W / 'preparation_v1/RESULT.json'
    assert pin(prep_path)['sha256'] == '24d4d3e68f38f402b49e9bccce05853384a3a6ddfdf929176527f63dc3693fa6'
    prep = json.loads(prep_path.read_text())
    config = json.loads(Path(prep['config']['path']).read_text())
    assert pin(prep['config']['path']) == prep['config']
    for name in ('training_result', 'objective', 'body_helper'):
        assert pin(config[name]['path']) == config[name]
    data = json.loads(Path(config['training_result']['path']).read_text())
    reviewer = load_module(D / 'selected_review_v1/REVIEW.py', 'root_selected_reviewer_initial_control')
    objective = load_module(config['objective']['path'], 'joint_objective_initial_control')
    body = load_module(config['body_helper']['path'], 'body_initial_control')
    model, saved_score = prep['initial']['model'], prep['initial']['score']
    exact = objective.evaluate(data, model, body)
    manual = reviewer.independent_score(data, model)
    checks = []
    def check(label, condition):
        assert condition, label
        checks.append(label)
    check('Exact frozen initial score replay', exact == saved_score)
    check('Complete independent residual denominator', len(manual['vector']) == len(exact['weighted_residuals']) == 413)
    deviations = [abs(a-b) for a,b in zip(manual['vector'], exact['weighted_residuals'])]
    check('All 413 residual weights/order agree', max(deviations) < 1e-12)
    check('Independent objective sum agrees', abs(manual['objective']-exact['objective']) < 1e-12)
    for group, rows in manual['energy_errors'].items():
        for name, value in rows.items():
            check('Energy sign/control ' + group + '/' + name,
                  abs(value-exact['energy_errors_kcal_mol'][group][name]) < 1e-12)
    for kind, rows in manual['force_checks'].items():
        for name, value in rows.items():
            target = exact['per_fixture'][kind][name]
            check('Force normalization/gate ' + kind + '/' + name,
                  abs(value['RMS']-target['rigid_RMS_kcal_mol_A']) < 1e-12 and
                  abs(value['maximum']-target['maximum_body_component_error_kcal_mol_A']) < 1e-12 and
                  value['passed'] == target['passed'])
    check('All ten energy controls included', sum(map(len, manual['energy_errors'].values())) == 10)
    check('All fifteen force controls included', sum(map(len, manual['force_checks'].values())) == 15)
    check('Initial failed W9 energy remains failed', abs(manual['energy_errors']['balanced_exchange']['aquo9_central']) > 1)
    check('Initial failed W9 La/Dy force gates remain failed', all(not manual['force_checks'][m]['aquo9_central']['passed'] for m in ('La', 'Dy')))
    check('No average hides initial development failure', manual['passed'] == exact['development_checks_passed'] == False)
    result = {
        'status': 'passed', 'checks_passed': len(checks), 'checks': checks,
        'implementation': pin(__file__), 'reviewer': pin(D / 'selected_review_v1/REVIEW.py'),
        'initial_fixture': pin(prep_path), 'training': config['training_result'],
        'objective': config['objective'], 'body_helper': config['body_helper'],
        'initial_objective': manual['objective'], 'max_weighted_residual_difference': max(deviations),
        'residual_count': 413, 'energy_control_count': 10, 'force_control_count': 15,
        'live_fit_result_read': False, 'new_native_calls': 0, 'numerical_threads': 1,
        'wall_seconds': time.monotonic()-start, 'CPU_seconds': time.process_time()-cpu,
    }
    out = W / 'selected_review_v1/source_review'
    out.mkdir(parents=True, exist_ok=True)
    dest = out / 'INITIAL_ARITHMETIC.json'
    assert not dest.exists()
    dest.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('status','checks_passed','initial_objective','max_weighted_residual_difference','wall_seconds','CPU_seconds')}))
    print(json.dumps(pin(dest)))

if __name__ == '__main__':
    main()
