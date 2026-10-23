#!/usr/bin/env python3
"""Artificial/rational cone correspondence only. Never instantiate a solver."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import sys
import time
from types import SimpleNamespace
import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('comparison', HERE/'COMPARE.py')
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
core = C.module(C.pin(HERE.parent/'CORE.py'))


def artificial():
    def block(n, den=4):
        # All deterministic rationals, independent of any saved scientific point.
        r = np.array([((i % 7)-3)/64 for i in range(n)])
        J = np.array([[((i+2*j) % 11-5)/32 for j in range(12)] for i in range(n)])
        return {'r': r, 'J': J, 'denominator': den}
    blocks = {name: block(n) for name, n in zip(C.BLOCKS, (413, 102, 12, 4, 2, 2))}
    forces = {}
    for kind in ('La', 'Dy', 'paired'):
        for name, n in zip(core.GEOMETRIES, (1, 1, 1, 8, 9)):
            b = block(3+6*n); b['weights'] = np.r_[np.full(3, 1/6), np.full(6*n, 1/(12*n))]
            forces[kind+'/'+name] = b
    energy = block(10); energy['names'] = ['control'+str(i) for i in range(10)]
    return {'blocks': blocks, 'forces': forces, 'energy': energy, 'q_current': 1.,
            'bounds_lower': np.r_[np.full(12, -.025), 0.], 'bounds_upper': np.r_[np.full(12, .025), 1.],
            'center': np.full(12, .5), 'width': np.ones(12), 'lower': np.zeros(12), 'upper': np.ones(12)}


def main(out):
    tick, cpu = time.monotonic(), time.process_time(); checks = []
    def check(name, ok):
        if not ok: raise AssertionError(name)
        checks.append(name)
    def rejects(name, call):
        try: call()
        except (ValueError, TypeError): checks.append(name); return
        raise AssertionError(name+' unexpectedly accepted')
    # A guard makes accidental backend construction in these controls fatal.
    def forbidden(*args, **kwargs): raise AssertionError('Solver construction forbidden in source controls')
    C.clarabel.DefaultSolver = forbidden
    p = artificial(); m = C.assemble(p)
    check('21_SOC_and856_linear', sum(r['kind']=='SOC' for r in m['layout']) == 21 and sum(r['kind']=='nonnegative' for r in m['layout']) == 856)
    check('1817_rows13_variables_zero_P', m['A'].shape == (1817,13) and m['P'].nnz == 0)
    differences = []
    for number, (step, t) in enumerate(((0, 1), (.01, .75), (-.02, .5), (.025, .25))):
        x = np.r_[[step*(-1 if j%2 else 1) for j in range(12)], t]
        r = C.correspondence(p, m, x, core); differences.append(r['maximum_arithmetic_difference'])
        check('all851_rational_probe_'+str(number), r['original_constraint_count'] == 851)
    # Independent normalization check: denominator4 means radius2, not4.
    row = m['layout'][0]; s = (m['b']-m['A']@np.r_[np.zeros(12),1])[row['start']:row['stop']]
    check('old_block_half_scale', np.array_equal(s[1:], p['blocks']['old_EF']['r']/2))
    for kind, target in (('energy', 1.), ('force', 5.)):
        for sign in (-1, 1):
            q = copy.deepcopy(p)
            b = q['energy'] if kind == 'energy' else next(iter(q['forces'].values()))
            b['r'][0] = sign*target; b['J'][0] = 0
            qm = C.assemble(q); values = qm['b']-qm['A']@np.r_[np.zeros(12),1]
            name = ('energy/'+('upper' if sign>0 else 'lower')+'/control0' if kind=='energy' else 'force/'+('upper' if sign>0 else 'lower')+'/La/water_minus/0')
            rr = next(r for r in qm['layout'] if r['name']==name)
            check(kind+'_exact_boundary_'+str(sign), values[rr['start']] == 0)
    bad = copy.deepcopy(p); bad['blocks']['old_EF']['denominator'] = 0
    rejects('nonpositive_normalization_rejected', lambda: C.assemble(bad))
    x = np.r_[np.zeros(12), .5]
    fractions = C.fraction_reports(p, x, 'Solved', core)
    check('original_four_fractions', [s['fraction'] for s in fractions] == [1., .5, .25, .125])
    check('q_not_t_interpolation', fractions[1]['certificate']['q_fraction'] == .625 and fractions[1]['certificate']['q_fraction'] != .75**2)
    check('q_selection_epsilon', all(s['q_selection_epsilon']==1e-10 for s in fractions))
    check('no_native_or_selection', all(s['native_evaluated'] is False for s in fractions))
    almost = C.fraction_reports(p, x, 'AlmostSolved', core)
    check('reduced_status_not_promoted', all(s['certificate']['solver_success'] is False and s['certificate']['failed_solver_status_remains_failed'] for s in almost))
    bound = copy.deepcopy(p); bound['energy']['r'][0] = 1.; bound['energy']['J'][0] = 0.
    yes = C.fraction_reports(bound, x, 'Solved', core); no = C.fraction_reports(bound, x, 'AlmostSolved', core)
    check('strict_failed_status_witness', yes[0]['eligible'] and not no[0]['eligible'])
    rejects('negative_t_not_squared_into_candidate', lambda: C.fraction_reports(p, np.r_[np.zeros(12),-.5], 'Solved', core))
    rejects('nonfinite_return_rejected', lambda: C.fraction_reports(p, np.r_[np.zeros(12),np.nan], 'Solved', core))
    tiny = {'A':sparse.csc_matrix([[-1.],[1.]]), 'b':np.array([0.,1.]), 'c':np.array([1.]),
            'layout':[{'name':'lower','kind':'nonnegative','start':0,'stop':1},{'name':'upper','kind':'nonnegative','start':1,'stop':2}]}
    diag = C.diagnostics(tiny, SimpleNamespace(x=[0.],s=[0.,1.],z=[1.,0.]))
    check('independent_dual_sign_and_gap', diag['primal_residual_inf']==diag['dual_residual_inf']==diag['signed_gap_t']==diag['complementarity_dot']==0)
    wrong = C.diagnostics(tiny, SimpleNamespace(x=[0.],s=[0.,1.],z=[-1.,0.]))
    check('wrong_dual_sign_visible', wrong['dual_residual_inf']==2 and wrong['dual_cone_margins']['lower']==-1)
    inp = {'campaign_driver':{'id':'driver'},'campaign_configuration':{'id':'config'},'campaign_fingerprint':'fp','expected_terminal_path':'literal'}
    cfg = {'terminal_campaign':{'path':'literal'}}
    fit = {'status':'bounded_component_fit_complete','mode':'run','cost':{},'implementation':inp['campaign_driver'],'configuration':inp['campaign_configuration'],'configuration_fingerprint':'fp'}
    C.terminal_gate(fit,cfg,inp); check('completed_terminal_accepted',True)
    C.terminal_gate(dict(fit,status='failed_component_fit'),cfg,inp); check('failed_terminal_preserved_and_accepted',True)
    for status in ('started','prepared_component_fit_zero_native','running'):
        rejects('nonterminal_rejected_'+status, lambda status=status:C.terminal_gate(dict(fit,status=status),cfg,inp))
    rejects('different_campaign_rejected', lambda:C.terminal_gate(dict(fit,configuration={'id':'other'}),cfg,inp))
    rejects('different_path_rejected', lambda:C.terminal_gate(fit,{'terminal_campaign':{'path':'other'}},inp))
    declared = json.loads((HERE/'CONFIG_UNBOUND.json').read_text())
    rejects('unbound_runtime_rejected', lambda:C.context(declared))
    settings=C.frozen_settings(C.pin(HERE/'SETTINGS.json'))
    check('all_settings_exact_and_bounded', settings.max_threads==1 and settings.time_limit==60 and settings.direct_solve_method=='qdldl')
    result={'status':'passed_artificial_conic_source_controls','implementation':C.pin(HERE/'COMPARE.py'),'test_source':C.pin(__file__),
            'checks':checks,'count':len(checks),'maximum_correspondence_difference':max(differences),'solver_constructions':0,'solver_calls':0,
            'native_calls':0,'scientific_vectors_read':0,'campaign_outcomes_read':False,'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}
    out=Path(out);out.mkdir(parents=True,exist_ok=False);C.save(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'checks':len(checks),'result':C.pin(out/'RESULT.json')}))


if __name__=='__main__': main(sys.argv[1])
