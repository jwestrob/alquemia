#!/usr/bin/env python3
"""Saved baseline and small algebra controls; no native or candidate evaluation."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def read(q):
    assert pin(q['path']) == q, q['path']
    return json.loads(Path(q['path']).read_text())


def module(q, name):
    assert pin(q['path']) == q
    s = importlib.util.spec_from_file_location(name, q['path'])
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
    return m


def save(p, value):
    Path(p).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    p = argparse.ArgumentParser(); p.add_argument('--config', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True); args = p.parse_args()
    assert all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'))
    args.out.mkdir(parents=True, exist_ok=False)
    tick, cpu = time.monotonic(), time.process_time(); checks = []
    def check(name, value):
        assert value, name
        checks.append(name)
    def close(a, b): return math.isclose(a, b, rel_tol=2e-14, abs_tol=2e-14)
    cfg = json.loads(args.config.read_text()); m = module(cfg['implementation'], 'joint_objective')
    inv = read(cfg['target_inventory']); data = read(cfg['development_data'])
    old = module(cfg['old_objective'], 'unchanged_objective'); body = module(cfg['body_helper'], 'unchanged_body')
    assert pin(cfg['declaration']['path']) == cfg['declaration'] and pin(cfg['specification']['path']) == cfg['specification']
    check('source helper pins equal target inventory', cfg['old_objective'] == inv['artifacts']['development_objective'] and cfg['body_helper'] == inv['artifacts']['body_helper'] and cfg['development_data'] == inv['artifacts']['development_data'])
    targets = m.targets_from_inventory(inv)
    response = {'pair_dipoles_Debye': {r['id']: r['selected17676_pair_mu_Debye'] for r in inv['groups']['pair_and_collective_dipoles']['pair_cycles']},
                'C_dipoles_Debye': {p: inv['groups']['pair_and_collective_dipoles']['parents'][p]['C']['selected17676_model_mu_Debye'] for p in m.PARENTS},
                'chi_states_a0_cubed': {r['quantity']: r['model_a0_cubed'] for r in targets['field_rows'] if r['quantity'] in m.STATES and r['step'] == 'h' and r['estimator'] == 'energy'}}
    original_response = copy.deepcopy(response); original_targets = copy.deepcopy(targets)
    result = m.evaluate(data, data['baseline_model'], response, targets, old, body)
    reference_old = old.evaluate(data, data['baseline_model'], body)
    check('old EF complete result returned unchanged', result['old_EF'] == reference_old)
    check('saved413 residuals exactly replayed', reference_old['weighted_residuals'] == data['baseline_score']['weighted_residuals'])
    check('original25 gates unchanged', result['old_EF']['energy_controls'] == data['baseline_score']['energy_controls'] and result['old_EF']['per_fixture'] == data['baseline_score']['per_fixture'])
    check('inputs not mutated', response == original_response and targets == original_targets)
    manual_pair = sum(sum(sum((response['pair_dipoles_Debye'][r['id']][j]-r['reference_Debye'][j])**2 for j in range(3)) for r in targets['pair_rows'] if r['parent'] == parent)/(3*(8 if parent.startswith('aquo8_') else 9)) for parent in m.PARENTS)/4
    manual_C = sum((response['C_dipoles_Debye'][p][j]-targets['C_Debye'][p][j])**2 for p in m.PARENTS for j in range(3))/12
    manual_chi = sum((response['chi_states_a0_cubed'][s]-targets['chi_a0_cubed'][s])**2 for s in m.STATES)/4
    check('saved baseline three response blocks match separate sums', all(close(result['block_objectives'][k], v) for k, v in [('pair', manual_pair), ('C', manual_C), ('chi', manual_chi)]))
    check('declared four equal block objective', close(result['objective'], (reference_old['objective']+manual_pair+manual_C+manual_chi)/4))
    check('531 weighted residual norm matches objective', len(result['weighted_residuals']) == 531 and close(sum(x*x for x in result['weighted_residuals']), result['objective']))
    check('exact selected four fixed field estimators', sum(r['used_in_objective'] for r in result['response']['field_rows']) == 4 and all(r['field_size_au'] == .001 and r['estimator'] == 'energy' and r['kind'] == 'state' for r in result['response']['field_rows'] if r['used_in_objective']))
    perfect = {'pair_dipoles_Debye': {r['id']: list(r['reference_Debye']) for r in targets['pair_rows']}, 'C_dipoles_Debye': copy.deepcopy(targets['C_Debye']), 'chi_states_a0_cubed': dict(targets['chi_a0_cubed'])}
    model_reference = {cid: {'energy': row['reference']['energy'], 'forces': row['reference']['forces']} for cid, row in data['rows'].items()}
    zero = m.evaluate(data, model_reference, perfect, targets, old, body)
    check('synthetic exact references give zero531 loss', zero['objective'] == 0 and all(x == 0 for x in zero['weighted_residuals']))
    for parent in (m.PARENTS[0], m.PARENTS[2]):
        q = copy.deepcopy(perfect)
        for row in targets['pair_rows']:
            if row['parent'] == parent: q['pair_dipoles_Debye'][row['id']][0] += 1
        check('equal parent weight despite8/9 pairs '+parent, close(m.response_blocks(q, targets)['objectives']['pair'], 1/12))
    q = copy.deepcopy(perfect); q['C_dipoles_Debye'][m.PARENTS[0]][0] += 1
    check('one C Cartesian error weight1/12', close(m.response_blocks(q, targets)['objectives']['C'], 1/12))
    q = copy.deepcopy(perfect); q['chi_states_a0_cubed'][m.STATES[0]] += 1
    check('one chi state error weight1/4', close(m.response_blocks(q, targets)['objectives']['chi'], .25))
    changed = copy.deepcopy(inv)
    for row in changed['groups']['conditional_susceptibility']['reported_rows']:
        if not (row['quantity'] in m.STATES and row['step'] == 'h' and row['estimator'] == 'energy'): row['reference_a0_cubed'] += 1000
    changed['groups']['trimer_energy_and_dipoles'] = {'excluded': 'not read'}
    for row in changed['groups']['pair_and_collective_dipoles']['parents'].values():
        row['P'] = row['B'] = {'excluded': 'not read'}
    changed_targets = m.targets_from_inventory(changed)
    check('derived field/P/B/E3 data cannot alter objective', m.response_blocks(response, changed_targets)['objectives'] == result['response']['objectives'])
    check('all other24 field rows retained descriptively', len(result['response']['field_rows']) == 28 and sum(not r['used_in_objective'] for r in result['response']['field_rows']) == 24)
    bad_cases = []
    q = copy.deepcopy(perfect); del q['pair_dipoles_Debye'][next(iter(q['pair_dipoles_Debye']))]; bad_cases.append(('missing pair', q))
    q = copy.deepcopy(perfect); del q['C_dipoles_Debye'][m.PARENTS[0]]; bad_cases.append(('missing C', q))
    q = copy.deepcopy(perfect); del q['chi_states_a0_cubed'][m.STATES[0]]; bad_cases.append(('missing chi', q))
    q = copy.deepcopy(perfect); q['chi_states_a0_cubed'][m.STATES[0]] = float('nan'); bad_cases.append(('nonfinite chi', q))
    q = copy.deepcopy(perfect); q['pair_dipoles_Debye'][next(iter(q['pair_dipoles_Debye']))] = [1, 2]; bad_cases.append(('bad vector', q))
    q = copy.deepcopy(perfect); q['B'] = {}; bad_cases.append(('undeclared additional response block', q))
    q = copy.deepcopy(perfect); q['chi_states_a0_cubed'][m.STATES[0]] = 1e300; bad_cases.append(('overflow block', q))
    for name, q in bad_cases:
        try: m.response_blocks(q, targets)
        except (ValueError, OverflowError): checks.append('reject '+name)
        else: raise AssertionError('Accepted '+name)
    save(args.out/'TARGETS.json', targets); save(args.out/'BASELINE_RESPONSE.json', response)
    save(args.out/'BASELINE_SCORE.json', result)
    receipt = {'status': 'passed_saved_objective_controls', 'implementation': cfg['implementation'], 'control_implementation': pin(__file__), 'configuration': pin(args.config), 'source_authorities': cfg,
               'targets': pin(args.out/'TARGETS.json'), 'baseline_response': pin(args.out/'BASELINE_RESPONSE.json'), 'baseline_score': pin(args.out/'BASELINE_SCORE.json'),
               'checks': checks, 'count': len(checks), 'saved_baseline_block_objectives': result['block_objectives'], 'saved_baseline_objective': result['objective'],
               'new_parameter_points_evaluated': 0, 'native_calls': 0, 'fit_performed': False,
               'cost': {'wall_seconds': time.monotonic()-tick, 'CPU_seconds': time.process_time()-cpu, 'threads': 1}}
    save(args.out/'RESULT.json', receipt)
    print(json.dumps({'result': pin(args.out/'RESULT.json'), 'checks': len(checks), 'baseline_blocks': result['block_objectives'], 'cost': receipt['cost']}))


if __name__ == '__main__':
    main()
