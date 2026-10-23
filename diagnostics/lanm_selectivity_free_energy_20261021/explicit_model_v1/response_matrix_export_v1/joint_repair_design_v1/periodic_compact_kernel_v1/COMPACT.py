"""Mathematical compact response prototype; no engine, geometry-file or fit API."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import time

BASE = Path(__file__).resolve().parent
TAIL = BASE.parent / 'periodic_tail_definition_v1/calculator_v1/TAIL.py'
TAIL_SHA = '875608169f73e618d55c1b1571294135d18818b14936d5885cc4260c78339136'
DECLARATION_SHA = '8c041a29b2631fdc8ab444eb50988d3fa4c43d851af929ccc3dcba3bace07263'


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def arithmetic_module():
    if pin(TAIL)['sha256'] != TAIL_SHA:
        raise ValueError('Outward arithmetic source changed')
    name = 'compact_kernel_pinned_tail_arithmetic'
    spec = importlib.util.spec_from_file_location(name, TAIL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def switch(r):
    """Return S,S',S''; accepts float or exact Fraction radii in Angstrom."""
    if r < 0:
        raise ValueError('Negative radius')
    z = r * 0
    if r <= 20:
        return z + 1, z, z
    if r >= 21:
        return z, z, z
    t = r - 20
    # Positive factored form of 1-10t^3+15t^4-6t^5.
    return ((1-t)**3 * (1+3*t+6*t*t),
            -30*t*t*(1-t)**2,
            -60*t*(1-t)*(1-2*t))


def radial(r, compact=True):
    """Pure scalar F,F',F'' prototype. Floating roundoff is not interval-certified."""
    if not math.isfinite(r) or r < 0:
        raise ValueError('Finite nonnegative radius required')
    s, sp, spp = switch(r) if compact else (1., 0., 0.)
    if s == sp == spp == 0:
        return 0., 0., 0.
    p = 1 + 2*r + 4*r*r/3
    f = 9*math.exp(10)/1849 * p*p * math.exp(-4*r)
    den = 3 + 6*r + 4*r*r
    num = 8*r + 16*r*r
    logp = -num/den
    logpp = -((8+32*r)*den-num*(6+8*r))/(den*den)
    fp, fpp = f*logp, f*(logp*logp+logpp)
    return f*s, fp*s+f*sp, fpp*s+2*fp*sp+f*spp


def norm(v):
    return math.sqrt(math.fsum(x*x for x in v))


def nearest(v, sides):
    if len(v) != 3 or len(sides) != 3 or any(not math.isfinite(x) for x in (*v, *sides)):
        raise ValueError('Three finite coordinates and cell sides required')
    if min(sides) < 44:
        raise ValueError('Outside compact prototype cell domain')
    return [x-L*math.floor(x/L+0.5) for x, L in zip(v, sides)]


def regularized_v(v, kappa):
    """Unchanged contact-regularized vector and Cartesian Jacobian."""
    r = norm(v)
    z = kappa*r**3
    if z < 1e-4:
        phi = kappa*math.fsum((-z)**n/math.factorial(n+1) for n in range(7))
        phip = 3*kappa*kappa*r*r*math.fsum(
            n*(-1)**n*z**(n-1)/math.factorial(n+1) for n in range(1, 7))
    else:
        s3 = -math.expm1(-z)
        phi = s3/r**3
        phip = 3*(z*math.exp(-z)-s3)/r**4
    value = [phi*x for x in v]
    jac = [[phi*int(i == j)+(phip/r*v[i]*v[j] if r else 0.)
            for j in range(3)] for i in range(3)]
    return value, jac


def response(vectors, *, a, b, c, d, kappa, alpha):
    """h, beta and their displacement Jacobians for artificial/prototype inputs.

    vectors are one compact nearest image per primary water. Both mutual RHS
    receive this same h; H has beta*I once on the primary metal. This routine
    does not solve, evaluate an energy, read coefficients, or certify stability.
    Metal coordinate derivatives are minus the sum of displacement derivatives.
    """
    vals = [a, b, c, d, kappa, alpha]
    if any(not math.isfinite(x) for x in vals) or kappa <= 0 or alpha <= 0:
        raise ValueError('Finite coefficients and positive alpha/kappa required')
    if len(vectors) > 3980 or any(len(v) != 3 or any(not math.isfinite(x) for x in v) for v in vectors):
        raise ValueError('Finite three-vectors within population cap required')
    radii = [norm(v) for v in vectors]
    ff = [radial(r) for r in radii]
    f = [x[0] for x in ff]
    df = [[fp*v[j]/r if r else 0. for j in range(3)]
          for v, r, (_, fp, _) in zip(vectors, radii, ff)]
    s = math.fsum(f)
    q = [math.fsum(f[j] for j in range(len(f)) if j != i) for i in range(len(f))]
    Q = math.fsum(x*x for x in f)
    collective = 2*math.fsum(f[i]*f[j] for i in range(len(f)) for j in range(i))
    eta, g = collective/(1+s*s), 2*s/(1+s)
    sigma = g*(c+d*eta)
    beta = math.expm1(sigma)/alpha
    h, dh, dbeta = [], [], []
    for i, v in enumerate(vectors):
        V, JV = regularized_v(v, kappa)
        A = a+b*q[i]/(1+q[i])
        h.append([3*f[i]*A*x for x in V])
        dh.append([[[3*A*(f[i]*JV[k][ell]+V[k]*df[i][ell]) if i == j
                     else 3*f[i]*b/(1+q[i])**2*V[k]*df[j][ell]
                    for ell in range(3)] for k in range(3)] for j in range(len(f))])
        deta = 2*q[i]/(1+s*s)-2*s*collective/(1+s*s)**2
        coeff = math.exp(sigma)/alpha*(2*(c+d*eta)/(1+s)**2+g*d*deta)
        dbeta.append([coeff*x for x in df[i]])
    return {'f': f, 'df': df, 's': s, 'Q': Q, 'q': q, 'eta': eta,
            'h': h, 'beta': beta, 'dh': dh, 'dbeta': dbeta}


def scalar_bounds():
    """Outward bounds, fixed uniform domain; no physical geometry evaluation."""
    tail = arithmetic_module()
    ar = tail.Arithmetic(precision=80)
    N, a, half, R = ar.integer(3980), ar.integer(44), ar.integer(22), ar.integer(20)
    rho = ar.exp(ar.integer(-176))
    # P(44+88j) has positive exact/interval coefficients.
    x, y = ar.integer(44), ar.integer(88)
    P = [ar.add(ar.add(ar.one, x), ar.div(ar.power(x, 2), ar.integer(3))),
         ar.add(y, ar.div(ar.mul(ar.integer(2), ar.mul(x, y)), ar.integer(3))),
         ar.div(ar.power(y, 2), ar.integer(3))]
    polynomial = tail.poly_mul(ar, [ar.integer(26), ar.integer(48), ar.integer(24)],
                               tail.poly_mul(ar, P, P))
    prefactor = ar.mul(N, ar.mul(tail.normalization(ar), ar.exp(ar.integer(-88))))
    shell_results = []
    for p in (0, 1):
        poly = polynomial if not p else tail.poly_mul(ar, polynomial, [half, a])
        if any(c.lo < 0 for c in poly):
            raise ValueError('Positive polynomial proof violated')
        series = ar.total(ar.mul(c, tail.geometric_moment(ar, rho, 0, j))
                          for j, c in enumerate(poly))
        shell_results.append(ar.mul(prefactor, series))
    M0, M1 = shell_results
    F20 = tail.radial_f(ar, R)
    nearest0 = ar.mul(N, F20)
    nearest1 = ar.mul(R, nearest0)
    B = ar.div(ar.integer(47), ar.integer(8))
    expressions = {
        'T0_weight_error': ('dimensionless', ar.add(nearest0, M0)),
        'T1_radial_derivative_error': ('angstrom^-1', ar.add(ar.mul(B, nearest0), ar.mul(ar.integer(4), M0))),
        'V0_radius_weighted_error': ('angstrom', ar.add(nearest1, M1)),
        'V1_unit_strain_derivative_error': ('dimensionless', ar.add(ar.mul(B, nearest1), ar.mul(ar.integer(4), M1))),
    }
    return {
        'domain': {'orthorhombic': True, 'minimum_side_angstrom': '44',
                   'maximum_primary_waters': 3980, 'switch_start_angstrom': '20',
                   'support_end_angstrom': '21', 'no_upper_cell_size_required': True},
        'arithmetic': {'precision': 80, 'Emin': -999999, 'Emax': 999999,
                       'outward_rounding': True, 'positive_geometric_series': True,
                       'fail_on_subnormal_underflow_overflow_clamping': True},
        'decomposition': {'F20': F20.json(), 'nearest_weight_bound': nearest0.json(),
                          'nearest_radius_weighted_bound': nearest1.json(),
                          'nonnearest_weight_bound': M0.json(),
                          'nonnearest_radius_weighted_bound': M1.json(),
                          'rho': rho.json(), 'polynomial_degrees': [6, 7],
                          'first_series_index': 0},
        'bounds': {key: {'units': units, 'positive_bound_expression_enclosure': val.json(),
                        'physical_error_enclosure': {'lower': '0', 'upper': str(val.hi)}}
                   for key, (units, val) in expressions.items()},
        'not_admitted': ['h/H coefficient-specific errors', 'self-consistent energy/force/stress',
                         'response stability', 'PME or native implementation', 'physical accuracy'],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    out = Path(args.out).resolve()
    if out.exists():
        raise ValueError('Immutable output already exists')
    if pin(BASE/'DECLARATION.md')['sha256'] != DECLARATION_SHA:
        raise ValueError('Prospective declaration changed')
    start_wall, start_cpu = time.perf_counter(), time.process_time()
    result = {'status': 'compact_kernel_scalar_bounds_failed',
              'sources': {'implementation': pin(__file__), 'declaration': pin(BASE/'DECLARATION.md'),
                          'outward_arithmetic': pin(TAIL)}, 'native_calls': 0,
              'physical_box_evaluations': 0, 'response_solves': 0}
    try:
        result['calculation'] = scalar_bounds()
        result['status'] = 'compact_kernel_scalar_bounds_complete'
    except Exception as e:
        result['error'] = repr(e)
        raise
    finally:
        result['cost'] = {'wall_seconds': time.perf_counter()-start_wall,
                          'process_CPU_seconds': time.process_time()-start_cpu,
                          'numerical_threads': 1}
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open('x') as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write('\n')


if __name__ == '__main__':
    main()
