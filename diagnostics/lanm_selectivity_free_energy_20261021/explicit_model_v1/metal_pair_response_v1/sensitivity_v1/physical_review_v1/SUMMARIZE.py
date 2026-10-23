#!/usr/bin/env python3
"""Finite saved-outcome interpretation; no native admission, fitting or calls."""
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[6]
BASE = Path('workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/metal_pair_response_v1/sensitivity_v1')
OUT = ROOT / BASE / 'physical_review_v1'


def pin(p):
    p = p.resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def norm(x):
    return math.sqrt(sum(v*v for v in x))


def compact_score(s):
    e = {f'{g}/{k}': {'error_kcal_mol': v['error_kcal_mol'], 'passed': v['passed']}
         for g, rows in s['energy_controls'].items() for k, v in rows.items()}
    f = {f'{g}/{k}': {a: v[a] for a in ('rigid_RMS_kcal_mol_A', 'maximum_body_component_error_kcal_mol_A', 'passed')}
         for g, rows in s['per_fixture'].items() for k, v in rows.items()}
    assert len(e) == 10 and len(f) == 15
    return {'objective': s['objective'], 'blocks': s['blocks'], 'energy': e, 'force': f,
            'energy_pass_count': sum(v['passed'] for v in e.values()),
            'force_pass_count': sum(v['passed'] for v in f.values()),
            'total_pass_count': sum(v['passed'] for v in (*e.values(), *f.values())),
            'all25_pass': s['development_checks_passed']}


def main():
    rp = ROOT / BASE / 'run_v1/execution/RESULT.json'
    mp = ROOT / BASE / 'preparation_v2/MANIFEST.json'
    r, m = json.loads(rp.read_text()), json.loads(mp.read_text())
    assert r['manifest'] == pin(mp)
    assert r['status'] == 'consumed_pair_response_sensitivity_complete'
    assert r['error'] is None and r['global_stop'] is False
    assert r['native_program_invocations'] == 236 and len(r['rows']) == 120
    assert all(v['status'] == 'complete' for v in r['rows'].values())
    assert r['internal_FD_energies'] == 0 and r['reused_gradient_programs'] == 4
    review_path = ROOT / BASE / 'independent_actual_review_v2/REVIEW.json'
    review_pin = pin(review_path)
    assert review_pin['sha256'] == 'd3c624ab02e1931905c6cb53ce7ccc582df0bf13187c0b20726a0d93bd81c379'
    baseline = compact_score(m['baseline_score'])
    result = {'status': 'saved_physical_interpretation_complete', 'implementation': pin(Path(__file__)),
              'application_result': pin(rp), 'manifest': pin(mp), 'baseline': baseline, 'profiles': {},
              'native_admission': review_pin,
              'admission_note': 'Independent saved review admitted all raw rows and arithmetic; only four proven tuple/list increments fields needed saved-review normalization. Executed source/results stay unchanged.',
              'new_native_calls': 0, 'fit_or_candidate_selection': False, 'full_shell_B_computed': False,
              'application_cost': r['cost'],
              'scheduler_cost_reported_by_root': {'job_id': '18413', 'elapsed_seconds': 48, 'allocated_CPUs': 4,
                 'allocated_CPU_seconds': 192, 'TotalCPU_seconds': 145.272, 'MaxRSS_K': 1493344}}
    for name, p in r['profiles'].items():
        s = compact_score(p['original25_controls'])
        s['increments'] = p['increments']
        s['objective_change'] = s['objective'] - baseline['objective']
        s['gained_checks'], s['lost_checks'] = [], []
        for category in ('energy', 'force'):
            for k, v in s[category].items():
                old = baseline[category][k]['passed']
                if old != v['passed']:
                    s['gained_checks' if v['passed'] else 'lost_checks'].append(f'{category}/{k}')
        s['trimer_cycles'] = {}
        for gid, t in p['trimer_cycles'].items():
            bmu, qmu, mu = t['baseline_mu3_Debye'], t['reference_mu3_Debye'], t['model_mu3_Debye']
            s['trimer_cycles'][gid] = {
                **{k: v for k, v in t.items() if k not in ('corrected_M_mu3_Debye', 'M_consistency_tolerance_Debye')},
                'E3_change_kcal_mol': t['model_E3_kcal_mol'] - t['baseline_E3_kcal_mol'],
                'baseline_E3_error_kcal_mol': t['baseline_E3_kcal_mol'] - t['reference_electronic_E3_kcal_mol'],
                'baseline_mu3_norm_Debye': norm(bmu), 'reference_mu3_norm_Debye': norm(qmu),
                'mu3_norm_change_Debye': norm(mu) - norm(bmu),
                'mu3_vector_error_norm_Debye': math.dist(mu, qmu),
                'baseline_mu3_vector_error_norm_Debye': math.dist(bmu, qmu),
                'mu3_magnitude_shortfall_percent': 100*(1-norm(mu)/norm(qmu))}
        t = list(s['trimer_cycles'].values())
        assert len(t) == 8
        s['trimer_summary'] = {
            'max_abs_E3_change_kcal_mol': max(abs(x['E3_change_kcal_mol']) for x in t),
            'E3_error_range_kcal_mol': [min(x['model_minus_reference_E3_kcal_mol'] for x in t), max(x['model_minus_reference_E3_kcal_mol'] for x in t)],
            'E3_error_improvement_count': sum(abs(x['model_minus_reference_E3_kcal_mol']) < abs(x['baseline_E3_error_kcal_mol']) for x in t),
            'mu3_vector_error_improvement_count': sum(x['mu3_vector_error_norm_Debye'] < x['baseline_mu3_vector_error_norm_Debye'] for x in t),
            'max_abs_mu3_norm_change_Debye': max(abs(x['mu3_norm_change_Debye']) for x in t),
            'mu3_magnitude_shortfall_percent_range': [min(x['mu3_magnitude_shortfall_percent'] for x in t), max(x['mu3_magnitude_shortfall_percent'] for x in t)]}
        s['paired_Dy_minus_La'] = {}
        for n in (8, 9):
            for where in ('near', 'far'):
                la = s['trimer_cycles'][f'aquo{n}_central_La__{where}']
                dy = s['trimer_cycles'][f'aquo{n}_central_Dy__{where}']
                pair = {}
                for label, ek, mk in [('baseline', 'baseline_E3_kcal_mol', 'baseline_mu3_norm_Debye'),
                                     ('model', 'model_E3_kcal_mol', 'model_mu3_norm_Debye'),
                                     ('reference', 'reference_electronic_E3_kcal_mol', 'reference_mu3_norm_Debye')]:
                    pair[label] = {'E3_kcal_mol': dy[ek]-la[ek], 'mu3_magnitude_Debye': dy[mk]-la[mk]}
                pair['mu3_magnitude_order_agrees'] = pair['model']['mu3_magnitude_Debye']*pair['reference']['mu3_magnitude_Debye'] > 0
                pair['E3_order_agrees'] = pair['model']['E3_kcal_mol']*pair['reference']['E3_kcal_mol'] > 0
                s['paired_Dy_minus_La'][f'W{n}_{where}'] = pair
        s['centers'] = {k: {a: b for a, b in v.items() if a != 'forces'} for k, v in p['centers'].items()}
        result['profiles'][name] = s
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / 'SUMMARY.json'
    assert not dest.exists(), 'Immutable output already exists'
    dest.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(pin(dest)))


if __name__ == '__main__':
    main()
