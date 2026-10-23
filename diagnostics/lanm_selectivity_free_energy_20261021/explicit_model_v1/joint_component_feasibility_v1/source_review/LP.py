"""Bounded tangent screening only; no energy evaluator or candidate selection."""
import json
import warnings

import numpy as np
from scipy.optimize import OptimizeWarning, linprog

NAMES = ('old_EF', 'pair_mu', 'C_mu', 'chi', 'P_energy', 'B_energy')
MODES = {'components_only': (), 'protect_EF': ('old_EF',),
         'protect_all_original': NAMES[:4]}


def solve(blocks, normalized):
    z = np.asarray(normalized, dtype=float)
    if (set(blocks) != set(NAMES) or z.ndim != 1 or not np.isfinite(z).all()
            or np.any(z < -1e-12) or np.any(z > 1+1e-12)):
        raise ValueError('Invalid blocks or normalized origin')
    z = np.clip(z, 0, 1)
    lower, upper = np.maximum(-.025, -z), np.minimum(.025, 1-z)
    arrays, gradients, norms = {}, {}, {}
    for name in NAMES:
        r = np.asarray(blocks[name]['residual'], dtype=float)
        j = np.asarray(blocks[name]['jacobian'], dtype=float)
        if (r.ndim != 1 or j.shape != (len(r), len(z))
                or not np.isfinite(r).all() or not np.isfinite(j).all()):
            raise ValueError('Invalid finite residual/Jacobian for '+name)
        arrays[name] = (r, j)
        norms[name] = float(np.linalg.norm(r))
        gradients[name] = (r @ j)/(r @ r) if norms[name] > 1e-12 else None
    degenerate = [n for n in NAMES if gradients[n] is None]
    if degenerate:
        return {'status': 'degenerate_block_norm', 'degenerate_blocks': degenerate,
                'block_norms': norms, 'modes': {}, 'native_evaluations': 0}
    modes = {}
    for name, protected in MODES.items():
        names = ('P_energy', 'B_energy') + protected
        A = np.array([np.r_[gradients[n], 1.0 if n in NAMES[4:] else 0.0]
                      for n in names])
        b = np.zeros(len(names))
        c = np.r_[np.zeros(len(z)), -1.0]
        bounds = list(zip(lower, upper)) + [(0, None)]
        # SciPy forwards the documented native HiGHS thread option unchanged.
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=OptimizeWarning,
                                    message='Unrecognized options detected.*threads.*')
            result = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method='highs',
                             options={'threads': 1, 'primal_feasibility_tolerance': 1e-8,
                                      'dual_feasibility_tolerance': 1e-8})
        if not result.success or not np.isfinite(result.x).all():
            raise ValueError('Tangent LP failed: '+str(result.message))
        x, dx = result.x, result.x[:-1]
        y = result.ineqlin.marginals
        lm, um = result.lower.marginals, result.upper.marginals
        dual = float(b @ y + lower @ lm[:-1] + upper @ um[:-1])
        # t has lower bound zero and no upper bound; its bound contribution is zero.
        primal = float(c @ x)
        gap = primal-dual
        stationarity = c - A.T @ y - lm - um
        feasibility = max(0., float(np.max(A @ x-b)), float(np.max(lower-dx)),
                          float(np.max(dx-upper)), float(-x[-1]))
        dual_sign = max(0., float(np.max(y)), float(np.max(-lm)),
                        float(np.max(um)), float(abs(um[-1])))
        complementarity = max(float(np.max(abs(y*(b-A@x)))),
                              float(np.max(abs(lm[:-1]*(dx-lower)))),
                              float(np.max(abs(um[:-1]*(upper-dx)))),
                              float(abs(lm[-1]*x[-1])))
        diagnostics = dict(primal_objective=primal, dual_objective=dual,
                           primal_minus_dual_gap=gap, primal_violation=feasibility,
                           dual_sign_violation=dual_sign,
                           stationarity_max=float(np.max(abs(stationarity))),
                           complementarity_max=complementarity)
        if (abs(gap) > 1e-7 or feasibility > 1e-7 or dual_sign > 1e-7
                or diagnostics['stationarity_max'] > 1e-7 or complementarity > 1e-7):
            raise ValueError('LP certificate failed: '+json.dumps(diagnostics))
        predictions = {}
        for n, (r, j) in arrays.items():
            predictions[n] = dict(origin_norm=norms[n],
                                  fractional_norm_directional_derivative=float(gradients[n]@dx),
                                  half_step_linearized_norm=float(np.linalg.norm(r+.5*j@dx)),
                                  full_step_linearized_norm=float(np.linalg.norm(r+j@dx)))
        modes[name] = dict(status='numerically_solved_tangent_problem',
                           solver_status=int(result.status), solver_message=result.message,
                           protected_blocks=list(protected), t=float(x[-1]),
                           normalized_direction=dx.tolist(), trust_lower=lower.tolist(),
                           trust_upper=upper.tolist(), constraint_names=list(names),
                           inequality_marginals=y.tolist(), lower_marginals=lm.tolist(),
                           upper_marginals=um.tolist(), certificate=diagnostics,
                           block_predictions=predictions,
                           interpretation='First-order aggregate screening; no finite-step or individual-gate guarantee.')
    return dict(status='complete_tangent_screen', block_norms=norms,
                fractional_norm_gradients={n:g.tolist() for n,g in gradients.items()},
                modes=modes, native_evaluations=0)


def selftest():
    def case(p, b, ef=(0., 0.)):
        return {n: {'residual': np.array([1.]),
                    'jacobian': np.array([p if n=='P_energy' else b if n=='B_energy'
                                          else ef if n=='old_EF' else (0., 0.)])}
                for n in NAMES}
    compatible = solve(case((-1., 0.), (-1., 0.)), [.5, .5])
    opposed = solve(case((1., 0.), (-1., 0.)), [.5, .5])
    protected = solve(case((-1., 0.), (-1., 0.), (1., 0.)), [.5, .5])
    boundary = solve(case((-1., 0.), (-1., 0.)), [1., .5])
    checks = {
        'compatible_analytic_optimum': abs(compatible['modes']['protect_all_original']['t']-.025)<1e-10,
        'opposed_zero_optimum': abs(opposed['modes']['components_only']['t'])<1e-10,
        'protected_conflict_zero': abs(protected['modes']['protect_EF']['t'])<1e-10,
        'domain_boundary_blocks_descent': abs(boundary['modes']['components_only']['t'])<1e-10,
        'full_linearized_norm': abs(compatible['modes']['components_only']['block_predictions']['P_energy']['full_step_linearized_norm']-.975)<1e-10,
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    return checks


if __name__ == '__main__':
    print(json.dumps(selftest(), indent=2))
