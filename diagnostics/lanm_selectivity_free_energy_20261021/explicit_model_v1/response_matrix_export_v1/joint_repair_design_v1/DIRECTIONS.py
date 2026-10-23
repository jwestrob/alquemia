"""Declared bounded linear proposals; no molecular evaluation or coefficient fit."""
import numpy as np
from scipy.optimize import lsq_linear


def directions(jacobian, residual, normalized):
    J, r, z = (np.asarray(x, dtype=float) for x in (jacobian, residual, normalized))
    if (J.shape != (len(r), len(z)) or not all(np.isfinite(x).all() for x in (J, r, z))
            or np.any(z < 0) or np.any(z > 1)):
        raise ValueError('Invalid finite linear problem or normalized box point')
    lower = np.maximum(-.025, -z)
    upper = np.minimum(.025, 1-z)
    full, _, rank, singular = np.linalg.lstsq(J, -r, rcond=1e-6)
    factor = 1.0
    for value, lo, hi in zip(full, lower, upper):
        if value > 0:
            factor = min(factor, hi/value)
        elif value < 0:
            factor = min(factor, lo/value)
    scaled = full*factor
    result = lsq_linear(J, -r, bounds=(lower, upper), method='bvls',
                        tol=1e-12, max_iter=1000, verbose=0)
    bounded = result.x
    if (not result.success or not np.isfinite(bounded).all()
            or np.any(bounded < lower-1e-12) or np.any(bounded > upper+1e-12)):
        raise ValueError('Declared bounded linear solver failed: '+str(result.message))
    # Raw gradient and active-set KKT residual are reported independently of
    # the library's optimality field. They are numerical diagnostics, not a
    # physical fitting tolerance or an alternative acceptance threshold.
    gradient = J.T@(J@bounded+r)
    kkt = gradient.copy()
    kkt[bounded <= lower+1e-12] = np.minimum(gradient[bounded <= lower+1e-12], 0)
    kkt[bounded >= upper-1e-12] = np.maximum(gradient[bounded >= upper-1e-12], 0)
    return scaled, bounded, {
        'rank': int(rank), 'singular_values': singular.tolist(),
        'unrestricted_GN_not_evaluated': full.tolist(), 'scaled_GN': scaled.tolist(),
        'bounded_LS': bounded.tolist(), 'bounded_solver': 'scipy.optimize.lsq_linear/bvls',
        'status': int(result.status), 'message': result.message, 'iterations': result.nit,
        'optimality': float(result.optimality), 'KKT_residual_max': float(np.max(abs(kkt))),
        'raw_gradient': gradient.tolist(), 'active_mask': result.active_mask.tolist(),
        'trust_lower': lower.tolist(), 'trust_upper': upper.tolist(),
    }
