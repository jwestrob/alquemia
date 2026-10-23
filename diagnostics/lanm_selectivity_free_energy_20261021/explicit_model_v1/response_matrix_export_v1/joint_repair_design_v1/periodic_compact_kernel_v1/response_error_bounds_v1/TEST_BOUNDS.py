"""Finite manufactured algebra controls; no actual amplification inputs."""
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import BOUNDS


def run():
    checks = []
    def check(name, passed):
        checks.append({'name': name, 'passed': bool(passed)})
        if not passed:
            raise AssertionError(name)
    def upper(result, key):
        return float(result['expressions'][key]['upper'])
    norms = dict(u_m='2', v_m='3', w_water_max='4', K_abs='1')
    spectral = dict(normalized_lower='0.5', alpha_max='2', alpha_o_max='1')
    zero = {k: '0' for k in ['T0', 'T1', 'V0', 'V1']}
    for direction in ['coordinate', 'strain']:
        out = BOUNDS.amplify(zero, direction, '1', '1.2', number_of_waters=3,
                             dual_norms=norms, spectral=spectral)
        check(direction+'_zero_difference', all(upper(out, k) == 0 for k in
              ['E_h', 'E_dh', 'E_beta', 'E_dbeta', 'E_normalized_H']))
        check(direction+'_zero_dipole_error', out['induced_dipoles']['status'] ==
              'conditional_perturbation_bound' and float(out['induced_dipoles']['E_u']['upper']) == 0)
    eps = {k: '0.001' for k in zero}
    missing = BOUNDS.amplify(eps, 'coordinate', '1', '1.2', number_of_waters=3)
    check('missing_dipoles_and_spectral_not_assessed', missing['explicit_gradient']['status'] ==
          missing['induced_dipoles']['status'] == 'not_assessed')
    bad = BOUNDS.amplify({k:'1' for k in zero}, 'coordinate', '1', '1.2',
                        number_of_waters=3, dual_norms=norms, spectral=spectral)
    check('unresolved_spectral_margin_inconclusive', bad['induced_dipoles']['status'] == 'inconclusive')
    for alpha_max, name in [('0.5', 'alpha_max_rejected')]:
        try:
            BOUNDS.amplify(eps, 'coordinate', '1', '1.2', number_of_waters=3,
                          dual_norms=norms, spectral={**spectral, 'alpha_max':alpha_max})
            check(name, False)
        except ValueError:
            check(name, True)

    # Independent literal environment/differential evaluation at manufactured
    # nonnegative weights. V,dV are bounded abstract vectors, not a physical box.
    def direct(f, df):
        s, ds = math.fsum(f), math.fsum(df)
        Q = math.fsum(x*x for x in f)
        dQ = 2*math.fsum(x*y for x,y in zip(f,df))
        eta = (s*s-Q)/(1+s*s)
        de = (2*s*(1+Q)*ds-(1+s*s)*dQ)/(1+s*s)**2
        g, dg = 2*s/(1+s), 2*ds/(1+s)**2
        a,b,c,d,alpha = .8,-.7,.1,-.15,1.2
        sigma, d_sigma = g*(c+d*eta), dg*(c+d*eta)+g*d*de
        beta, dbeta = math.expm1(sigma)/alpha, math.exp(sigma)*d_sigma/alpha
        h, dh = [], []
        for i in range(3):
            q = math.fsum(f[j] for j in range(3) if j != i)
            dq = math.fsum(df[j] for j in range(3) if j != i)
            V = [.3/(i+1), -.2, .1]
            dV = [.05, .03/(i+1), -.04]
            A = a+b*q/(1+q)
            h.append([3*f[i]*A*x for x in V])
            dh.append([3*(df[i]*A*V[j]+f[i]*A*dV[j]+f[i]*b*dq/(1+q)**2*V[j])
                       for j in range(3)])
        return h, dh, beta, dbeta
    fc, dfc = [.1,.2,.3], [-.02,.03,.01]
    delta, ddelta = [.01,.005,.02], [.001,-.002,.0005]
    for name, dd in [('general_attenuation',ddelta), ('retained_chain_only',[0.,0.,0.])]:
        full = direct([x+y for x,y in zip(fc,delta)], [x+y for x,y in zip(dfc,dd)])
        compact = direct(fc,dfc)
        T = math.fsum(delta)
        D = math.fsum(abs(x) for x in dd)
        result = BOUNDS.amplify({'T0':str(T),'T1':str(D),'V0':str(T),'V1':str(D)},
                'coordinate','1','1.2',number_of_waters=3,
                retained={'s':str(math.fsum(fc)), 'D':str(math.fsum(abs(x) for x in dfc)),
                          'W':str(math.fsum(fc))}, dual_norms=norms)
        errors = [math.fsum(math.sqrt(math.fsum((a-b)**2 for a,b in zip(x,y)))
                            for x,y in zip(full[k],compact[k])) for k in [0,1]]
        errors += [abs(full[k]-compact[k]) for k in [2,3]]
        for key, error in zip(['E_h','E_dh','E_beta','E_dbeta'],errors):
            check(name+'_'+key,error <= upper(result,key))
        # Fixed-dual triangle bound independently applied to literal differences.
        direct_error = .5*(6*errors[3]+4*errors[1])
        check(name+'_explicit_gradient',direct_error <=
              float(result['explicit_gradient']['error_expression']['upper']))
    return checks


if __name__ == '__main__':
    out=Path(sys.argv[1]).resolve()
    if out.exists(): raise ValueError('Immutable output exists')
    wall,cpu=time.perf_counter(),time.process_time()
    result={'status':'failed','native_calls':0,'physical_box_evaluations':0,
            'physical_response_solves':0,'actual_bound_amplifications':0,
            'source_sha256':hashlib.sha256(Path(BOUNDS.__file__).read_bytes()).hexdigest(),
            'test_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    try:
        result['checks']=run();result['check_count']=len(result['checks']);result['status']='passed'
    except Exception as e:
        result['error']=repr(e);raise
    finally:
        result['cost']={'wall_seconds':time.perf_counter()-wall,
                        'process_CPU_seconds':time.process_time()-cpu,'numerical_threads':1}
        out.parent.mkdir(parents=True,exist_ok=True)
        with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
