#!/usr/bin/env python3
"""Declared saved review: raw16 site files -> four baseline Jacobian rows.

Recompute both steps, native units, print bounds, ranks and all domain envelopes.
No native execution, fitted coefficients or repeated model19 admission.
"""
import datetime
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np

ROOT = Path.cwd()
BRANCH = 'lanm_selectivity_free_energy_20261021/explicit_model_v1/conditional_susceptibility_v1/comparison/saved_stiffness_jacobian_v1'
STATES = ['La_dimer', 'La_trimer', 'Dy_dimer', 'Dy_trimer']
TERMS = {'La': [(-1, 'La_dimer'), (1, 'La_trimer')],
         'Dy': [(-1, 'Dy_dimer'), (1, 'Dy_trimer')],
         'Dy_minus_La': [(1, 'La_dimer'), (-1, 'La_trimer'), (-1, 'Dy_dimer'), (1, 'Dy_trimer')]}

def pin(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())

def raw(p):
    assert pin(p['path']) == p
    return Path(p['path']).read_text()

def read(p):
    return json.loads(raw(p))

def main(out):
    tick, cpu = time.monotonic(), time.process_time()
    target = pin(ROOT/'workspaces'/BRANCH/'actual_v1/RESULT.json')
    actual = read(target)
    assert actual['status'] == 'saved_oxygen_stiffness_jacobian_complete'
    cfg = read(actual['configuration'])
    for key in ('implementation', 'declaration', 'design', 'model_review'):
        raw(cfg[key])
    m, native, comp = [read(cfg[k]) for k in ('model_manifest', 'model_result', 'comparison_result')]
    units = m['conversion']
    bohr, debye = map(Decimal, map(str, (units['bohr_A'], units['debye_per_eA'])))
    prefactor = float(Decimal(str(units['electric']))*bohr**2/Decimal(str(units['hartree_kcal_mol'])))
    records, diffs = {}, []
    def same(x, y):
        a, b = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        assert a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
        d = float(np.max(np.abs(a-b))) if a.size else 0.
        diffs.append(d)
        assert d <= 1e-10, d
    for state in STATES:
        row = m['rows'][state]
        xyz = [line.split() for line in raw(row['input']['xyz']).splitlines()[1:] if line.strip()]
        coords = np.array([[float(t) for t in line[2:5]] for line in xyz])
        assert [int(line[5]) for line in xyz] == row['types']
        alphas = {}
        for source in ('parameters', 'key'):
            for line in raw(row['input'][source]).splitlines():
                t = line.split('#', 1)[0].split()
                if t and t[0].lower() == 'polarize': alphas[int(t[1])] = float(t[2])
        assert alphas[36] == .837
        oxy = [i for i, t in enumerate(row['types']) if t == 36]
        distances = [float(np.linalg.norm(coords[i]-coords[0])) for i in oxy]
        normals = [(coords[i]-coords[0])/r for i, r in zip(oxy, distances)]
        f = [(1+2*r+(2*r)**2/3)**2/(1+5+25/3)**2*math.exp(-4*(r-2.5)) for r in distances]
        eta = [math.fsum(f[j] for j in range(len(f)) if i != j) for i in range(len(f))]
        records[state] = {}
        for tag, step in [('h', '.001'), ('half', '.0005')]:
            pm = []
            for sign in ('plus', 'minus'):
                rec = native['rows'][state+'__'+sign+'_'+tag]
                assert rec['status'] == 'complete'
                lines = [v.split() for v in raw(rec['moments']['induced_output']).splitlines() if v.strip()]
                assert int(lines[0][0]) == row['atoms'] and len(lines) == row['atoms']+1
                assert [int(v[0]) for v in lines[1:]] == list(range(1, row['atoms']+1))
                assert [v[1] for v in lines[1:]] == row['native_names']
                pm.append([[Decimal(s) for s in v[2:]] for v in lines[1:]])
            with localcontext() as context:
                context.prec = 40
                denom = 2*Decimal(step)*debye*bohr
                x = np.array([[float((a-b)/denom) for a,b in zip(u,v)] for u,v in zip(*pm)])
                eps = float(Decimal('.0000005')/(Decimal(step)*debye*bohr))
            projections = np.array([n@x[i] for n, i in zip(normals, oxy)])
            pbound = np.array([sum(abs(n))*eps for n in normals])
            squares = projections**2
            square_bound = (np.abs(projections)+pbound)**2-squares
            weight = prefactor/alphas[36]
            db = -weight*np.dot(f, squares)
            dc = -weight*np.dot(np.multiply(f, eta), squares)
            eb = weight*np.dot(f, square_bound)
            ec = weight*np.dot(np.multiply(f, eta), square_bound)
            J, err = np.zeros(6), np.zeros(6)
            j = 1 if state.startswith('La') else 4
            J[j:j+2], err[j:j+2] = [db, dc], [eb, ec]
            sigma, se = .2*weight*sum(squares), .2*weight*sum(square_bound)
            saved = actual['states'][state]['steps'][tag]
            same(J, saved['J_a0_cubed_per_parameter']); same(err, saved['J_print_bound'])
            same(x, saved['x_site_au']); same(sigma, saved['sigma_domain_first_order_abs_upper_bound'])
            same(se, saved['sigma_domain_upper_bound_print_allowance'])
            same(sum(x[:,2]), saved['chi_site_sum_a0_cubed'])
            records[state][tag] = dict(J=J.tolist(), print_bound=err.tolist(), sigma=float(sigma), sigma_print=float(se))
    for name, terms in TERMS.items():
        records[name] = {}
        for tag in ('h', 'half'):
            J = sum(c*np.array(records[s][tag]['J']) for c,s in terms)
            err = sum(abs(c)*np.array(records[s][tag]['print_bound']) for c,s in terms)
            sigma = sum(abs(c)*records[s][tag]['sigma'] for c,s in terms)
            se = sum(abs(c)*records[s][tag]['sigma_print'] for c,s in terms)
            same(J, actual['derived'][name]['steps'][tag]['J_a0_cubed_per_parameter'])
            records[name][tag] = dict(J=J.tolist(), print_bound=err.tolist(), sigma=sigma, sigma_print=se)
    rank = {}
    for tag in ('h', 'half'):
        full = np.array([records[s][tag]['J'] for s in STATES])
        rank[tag] = {}
        for family, columns in [('six', [0,1,2,3,4,5]), ('nested_four', [0,1,3,4])]:
            sv = np.linalg.svd(full[:,columns], compute_uv=False)
            saved = actual['rank'][tag][family]
            same(sv[:len(saved['nonzero_candidate_singular_values'])], saved['nonzero_candidate_singular_values'])
            n = int(sum(sv > 1e-6*sv[0]))
            assert n == saved['rank_at_relative_1e_6']
            rank[tag][family] = dict(rank=n, singular_values=sv.tolist())
    envelope_rows = []
    lookup = {(r['quantity'],r['step'],r['estimator'],r['family']):r for r in actual['linear_envelopes']}
    assert len(lookup) == len(actual['linear_envelopes']) == 56 and len(comp['table']) == 28
    for r in comp['table']:
        s = records[r['quantity']][r['step']]
        for family, columns in [('six', [1,2,4,5]), ('nested_four', [1,4])]:
            box = .2*sum(abs(s['J'][j]) for j in columns)
            box_print = .2*sum(abs(s['J'][j])+s['print_bound'][j] for j in columns)
            upper, printed = min(box,s['sigma']), min(box_print,s['sigma']+s['sigma_print'])
            need = -r['model_minus_reference_a0_cubed']
            outside = abs(need) > printed+r['difference_print_bound_a0_cubed']
            old = lookup[(r['quantity'],r['step'],r['estimator'],family)]
            same(upper,old['conservative_first_order_abs_upper_bound']);same(printed,old['upper_bound_including_site_print_allowance'])
            assert outside == old['required_shift_exceeds_first_order_envelope']
            envelope_rows.append(dict(quantity=r['quantity'],step=r['step'],estimator=r['estimator'],family=family,needed=need,upper=upper,outside=outside))
    state_checks = [x for x in envelope_rows if x['quantity'] in STATES]
    assert len(state_checks) == 32 and all(x['outside'] for x in state_checks)
    paired6 = [x for x in envelope_rows if x['quantity']=='Dy_minus_La' and x['family']=='six']
    assert len(paired6) == 4 and not any(x['outside'] for x in paired6)
    result = dict(status='passed_independent_saved_stiffness_arithmetic_review', review_source=pin(__file__), input_result=target,
                  configuration=actual['configuration'], model_actual_review=cfg['model_review'], raw_files_replayed=16,
                  states_and_derived=records, independent_SVD=rank, envelopes=envelope_rows,
                  maximum_arithmetic_difference=max(diffs), state_envelopes_exceeded=32, state_envelopes_declared=32,
                  paired_six_inside_bounds=4, paired_six_declared=4,
                  conclusion='All four state offsets exceed both-family first-order bounds at both steps/estimators. The six-family paired contrast alone lies inside a conservative bound; this neither meets the state constraints nor proves simultaneous feasibility. Stiffness rank is4 (nested2), so local amplitude leverage, not curvature-column rank loss, is the observed limitation. No finite-parameter impossibility or unique physics conclusion.',
                  actual_calculation_cost=actual['cost'], new_native_calls=0, full_model_admission_repeated=False,
                  review_cost=dict(wall_seconds=time.monotonic()-tick, CPU_seconds=time.process_time()-cpu, threads=1))
    out.mkdir(parents=True, exist_ok=False)
    path=out/'REVIEW.json';path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'], review=pin(path), maximum_difference=max(diffs))))

if __name__ == '__main__':
    main(ROOT/'workspaces'/BRANCH/'independent_actual_review_v1')
