"""Fixed supernodal LL and sparse sampled residuals; no engine or data loader."""
import warnings
import numpy as np
from scipy.sparse import csc_array
from sksparse import cholmod

SETTINGS = dict(lower=True, order='amd', sym_kind='sym', supernodal_mode='supernodal')
EPS = np.finfo(float).eps
TINY = np.finfo(float).tiny


def validate_lower(a):
    if not isinstance(a, csc_array) or a.dtype != np.dtype('float64'):
        raise ValueError('canonical float64 csc_array required')
    if a.ndim != 2 or a.shape[0] != a.shape[1] or a.shape[0] == 0:
        raise ValueError('nonempty square matrix required')
    if not a.has_canonical_format or not np.all(np.isfinite(a.data)):
        raise ValueError('noncanonical or nonfinite lower matrix')
    for j in range(a.shape[0]):
        rows = a.indices[a.indptr[j]:a.indptr[j+1]]
        if len(rows) == 0 or rows[0] != j:
            raise ValueError('explicit diagonal and lower-only structure required')
    # One triangle alone cannot establish symmetry of an upstream full operator.
    return a


def infnorm(v):
    return float(np.max(np.abs(v)))


def matrix_norm(a):
    z = abs(a)
    rows = np.asarray(z.sum(axis=1)).ravel()
    cols = np.asarray(z.sum(axis=0)).ravel()
    return float(np.max(rows + cols - np.abs(a.diagonal())))


def action(a, v):
    return np.asarray(a @ v + a.T @ v - a.diagonal()*v)


def metadata(f):
    return dict(N=int(f.N), is_numeric=bool(f.is_numeric), is_ll=bool(f.is_ll),
                is_super=bool(f.is_super), is_lower=bool(f.is_lower), order=f.order,
                raw_minor=None, raw_status=None,
                unavailable_reason='raw C minor/status are not exposed by this public wrapper')


def factor_once(a, on_phase=None):
    """Exactly one symbolic constructor and one numeric call, without fallback."""
    validate_lower(a)
    if on_phase is not None:
        on_phase('symbolic_analyses')
    f = cholmod.CholeskyFactor(a, **SETTINGS)
    error = None
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter('always')
        try:
            if on_phase is not None:
                on_phase('numeric_factorizations')
            f.factorize(a, ldl=False, beta=0.)
        except Exception as exc:
            error = dict(type=type(exc).__name__, message=str(exc),
                         attributes={k: repr(v) for k, v in vars(exc).items()})
    warning_rows = [dict(category=type(w.message).__name__, message=str(w.message)) for w in seen]
    return f, error, warning_rows


def ten_probes(n, metal_coordinates=(0, 1, 2)):
    """Fixed vectors; n<3 tiny fixtures use cyclic duplicate Cartesian slots."""
    if n < 1:
        raise ValueError('positive dimension required')
    vectors = []
    axis = np.arange(n) % min(3, n)
    for k in range(3):
        vectors.append((f'uniform_axis_{k}', (axis == k % min(3, n)).astype(float)))
    for k, coord in enumerate(metal_coordinates):
        if n >= 3 and not 0 <= coord < n:
            raise ValueError('metal coordinate outside active map')
        v = np.zeros(n); v[coord % n] = 1.
        vectors.append((f'metal_basis_{k}', v))
    rng = np.random.Generator(np.random.PCG64(20261022))
    for k in range(4):
        v = 2.*rng.integers(0, 2, size=n).astype(float) - 1.
        vectors.append((f'rademacher_{k}', v))
    if len(vectors) != 10:
        raise ValueError('exact ten-probe inventory required')
    return vectors


def factor_action_error(a, l, p, v, knorm):
    got = np.asarray(l @ (l.T @ v[p]))
    wanted = action(a, v)[p]
    err = infnorm(got-wanted)/(max(1., knorm)*infnorm(v))
    return err, bool(np.all(np.isfinite(got)) and np.all(np.isfinite(wanted)))


def audit_factor(f, a, probes, on_solve=None):
    """One copied sparse L; no dense matrix or matrix product L L-transpose."""
    validate_lower(a)
    n = a.shape[0]
    before = metadata(f)
    p = np.asarray(f.get_perm()).copy()
    if p.shape != (n,) or not np.array_equal(np.sort(p), np.arange(n)):
        raise ValueError('invalid permutation')
    l = f.get_factor(kind='LL', lower=True)
    after = metadata(f)
    knorm = matrix_norm(a)
    bound = 100*n*EPS
    positive = bool(l.shape == a.shape and np.all(np.isfinite(l.data))
                    and np.all(np.isfinite(l.diagonal())) and np.all(l.diagonal() > 0))
    rows = []
    for name, v in probes:
        if v.shape != (n,) or not np.all(np.isfinite(v)) or infnorm(v) == 0:
            raise ValueError('invalid probe')
        ferr, finite_action = factor_action_error(a, l, p, v, knorm)
        if on_solve is not None:
            on_solve()
        x = np.asarray(f.solve(v, system='A'))
        residual = action(a, x)-v
        serr = infnorm(residual)/max(TINY, knorm*infnorm(x)+infnorm(v))
        finite = bool(finite_action and np.all(np.isfinite(x))
                      and np.all(np.isfinite(residual)) and np.isfinite(ferr) and np.isfinite(serr))
        rows.append(dict(name=name, factor_action_error=ferr, solve_backward_error=serr,
                         finite=finite, passed=finite and ferr <= bound and serr <= bound))
    flags = bool(before['is_numeric'] and before['is_ll'] and before['is_super']
                 and before['is_lower'] and before == after)
    rec = dict(before_copy=before, after_copy=after, permutation=p.tolist(),
               sparse_factor_shape=list(l.shape), sparse_factor_nnz=int(l.nnz),
               finite_positive_diagonal=positive, original_factor_flags_preserved=flags,
               K_infinity_norm=knorm, arithmetic_bound=bound, probes=rows,
               maximum_factor_action_error=max(r['factor_action_error'] for r in rows),
               maximum_solve_backward_error=max(r['solve_backward_error'] for r in rows),
               passed=positive and flags and all(r['passed'] for r in rows),
               limitation='Sampled action and solve backward checks; not a whole factor-error norm or spectral forward-error certificate.')
    return rec, l, p
