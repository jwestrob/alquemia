"""Tiny artificial controls; no physical inputs, bound run, or response solve."""
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path
import sys
import time

import COMPACT as C


def run():
    checks = []
    def check(name, condition, **details):
        checks.append({'name': name, 'passed': bool(condition), **details})
        if not condition:
            raise AssertionError(name)
    def close(name, actual, expected, floor=1e-12):
        error = abs(actual-expected)
        allowance = 2e-7*max(abs(actual), abs(expected))+floor
        check(name, error <= allowance, absolute_error=error, allowance=allowance)

    for radius, expected in [(F(0), (1, 0, 0)), (F(20), (1, 0, 0)),
                             (F(21), (0, 0, 0)), (F(23), (0, 0, 0))]:
        check('exact_switch_boundary_'+str(radius), C.switch(radius) == expected)
    # Rational identities supplement the global factored sign proof.
    for t in [F(1, 10), F(1, 4), F(1, 2), F(3, 4), F(9, 10)]:
        s, sp, spp = C.switch(20+t)
        check('quintic_sign_and_derivative_'+str(t),
              s == 1-10*t**3+15*t**4-6*t**5 and 0 < s < 1
              and sp == -30*t*t*(1-t)**2 and -F(15, 8) <= sp < 0
              and spp == -60*t+180*t*t-120*t**3)
    check('slope_max_attained', C.switch(F(41, 2))[1] == -F(15, 8))
    for r in [F(0), F(1, 2), F(1), F(20), F(22), F(66)]:
        ratio = 8*r*(1+2*r)/(3+6*r+4*r*r)
        check('analytic_F_monotonic_'+str(r), 0 <= ratio <= 4)
        if r >= 20:
            check('analytic_rF_monotonic_'+str(r), r*ratio > 1)
    check('zero_kernel_matches_original_contact', C.radial(0) == C.radial(0, compact=False))
    check('C2_inner_join', C.radial(20) == C.radial(20, compact=False))
    check('C2_outer_join', C.radial(21) == (0., 0., 0.))
    h = 1e-5
    for r in (2.5, 20.25, 20.5, 20.75):
        f, fp, fpp = C.radial(r)
        minus, plus = C.radial(r-h), C.radial(r+h)
        close('radial_first_'+str(r), (plus[0]-minus[0])/(2*h), fp, floor=0.)
        close('radial_second_'+str(r), (plus[1]-minus[1])/(2*h), fpp, floor=0.)
    f, fp, _ = C.radial(20.5, compact=False)
    s, sp, _ = C.switch(20.5)
    correct = C.radial(20.5)[1]
    check('switch_derivative_is_material', abs(correct-fp*s) > abs(correct)*0.1)
    close('explicit_switch_product_rule', correct, fp*s+f*sp, floor=0.)

    # Coordinate and homogeneous-strain derivatives at a transition point.
    v = [20.2, 2.0, 0.4]
    r = C.norm(v)
    grad = [C.radial(r)[1]*x/r for x in v]
    for k in range(3):
        vm, vp = v.copy(), v.copy()
        vm[k] -= h; vp[k] += h
        close('transition_cartesian_'+str(k),
              (C.radial(C.norm(vp))[0]-C.radial(C.norm(vm))[0])/(2*h), grad[k], floor=0.)
    strain = [[0.2, 0.1, 0.0], [0.1, -0.15, 0.05], [0.0, 0.05, 0.1]]
    tangent = [math.fsum(strain[i][j]*v[j] for j in range(3)) for i in range(3)]
    vm = [v[i]-h*tangent[i] for i in range(3)]
    vp = [v[i]+h*tangent[i] for i in range(3)]
    close('transition_strain_virial',
          (C.radial(C.norm(vp))[0]-C.radial(C.norm(vm))[0])/(2*h),
          math.fsum(grad[i]*tangent[i] for i in range(3)), floor=0.)

    # Finite image enumeration is exact for these artificial compact cases:
    # omitted shells start at >=66 Angstrom, outside support.
    images = list(itertools.product(range(-1, 2), repeat=3))
    for index, (sides, raw) in enumerate([
            ([44., 44., 44.], [20.5, 0., 0.]),
            ([44., 50., 60.], [-20.5, 0., 0.]),
            ([44., 44., 44.], [22., 0., 0.]),
            ([44., 50., 60.], [65.5, -50., 60.]),
            ([44., 50., 60.], [0., 0., 0.])]):
        v = C.nearest(raw, sides)
        terms = [C.radial(C.norm([v[j]+n[j]*sides[j] for j in range(3)]))[0]
                 for n in images]
        nearest = C.radial(C.norm(v))[0]
        check('compact_lattice_identity_'+str(index), math.fsum(terms) == nearest
              and sum(x != 0 for x in terms) <= 1)
        shifted = [raw[j]+[3, -2, 1][j]*sides[j] for j in range(3)]
        check('reimage_invariance_'+str(index), C.nearest(shifted, sides) == v)
    for k in (1, 2, 3):
        actual = sum(max(abs(x) for x in n) == k
                     for n in itertools.product(range(-k, k+1), repeat=3))
        check('exact_shell_count_'+str(k), actual == 24*k*k+2)
    try:
        C.nearest([0., 0., 0.], [43.9, 44., 44.])
        check('invalid_domain_rejected', False)
    except ValueError:
        check('invalid_domain_rejected', True)

    # Artificial amplitude values only; no fitted or selected coefficients.
    params = {'a': .1, 'b': -.07, 'c': .08, 'd': -.03, 'kappa': .2, 'alpha': 1.3}
    vectors = [[2.4, .2, -.1], [-.3, 2.6, .4], [.2, -.4, 2.8]]
    base = C.response(vectors, **params)
    max_h_error = max_beta_error = 0.
    for j in range(3):
        for ell in range(3):
            vm, vp = [v.copy() for v in vectors], [v.copy() for v in vectors]
            vm[j][ell] -= h; vp[j][ell] += h
            rm, rp = C.response(vm, **params), C.response(vp, **params)
            db = (rp['beta']-rm['beta'])/(2*h)
            max_beta_error = max(max_beta_error, abs(db-base['dbeta'][j][ell]))
            close('environment_beta_'+str((j, ell)), db, base['dbeta'][j][ell])
            for i in range(3):
                for k in range(3):
                    dh = (rp['h'][i][k]-rm['h'][i][k])/(2*h)
                    max_h_error = max(max_h_error, abs(dh-base['dh'][i][j][k][ell]))
                    close('environment_h_'+str((i, j, k, ell)), dh, base['dh'][i][j][k][ell])
    tail_response = C.response([[21., 0., 0.], [0., 22., 0.]], **params)
    check('compact_all_environment_terms_zero', tail_response['s'] == tail_response['Q'] ==
          tail_response['beta'] == 0. and all(x == 0. for row in tail_response['h'] for x in row)
          and all(x == 0. for row in tail_response['dbeta'] for x in row))
    # No dependency instantiation/solver or real-box access occurs here.
    return {'checks': checks, 'check_count': len(checks),
            'max_artificial_h_derivative_error': max_h_error,
            'max_artificial_beta_derivative_error': max_beta_error}


if __name__ == '__main__':
    out = Path(sys.argv[1]).resolve()
    if out.exists():
        raise ValueError('Immutable control output already exists')
    wall, cpu = time.perf_counter(), time.process_time()
    result = {'status': 'failed', 'source': C.pin(C.__file__), 'controls': C.pin(__file__),
              'declaration': C.pin(C.BASE/'DECLARATION.md'), 'native_calls': 0,
              'physical_box_evaluations': 0, 'response_solves': 0}
    try:
        result.update(run())
        result['status'] = 'passed'
    except Exception as e:
        result['error'] = repr(e)
        raise
    finally:
        result['cost'] = {'wall_seconds': time.perf_counter()-wall,
                          'process_CPU_seconds': time.process_time()-cpu,
                          'numerical_threads': 1}
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open('x') as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write('\n')
