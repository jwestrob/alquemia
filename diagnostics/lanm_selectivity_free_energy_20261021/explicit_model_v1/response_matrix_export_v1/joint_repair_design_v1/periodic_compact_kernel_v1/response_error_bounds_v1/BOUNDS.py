"""Pure one-metal amplification contract; no geometry, fit, native or solve API."""
from pathlib import Path
import importlib.util
import hashlib
import sys

BASE = Path(__file__).resolve().parent
COMPACT_PATH = BASE.parent/'COMPACT.py'
COMPACT_SHA = 'f99522d8f967fe8a9761df7c588ebbd07fe43409d889784e761c0232069befae'


def arithmetic():
    if hashlib.sha256(COMPACT_PATH.read_bytes()).hexdigest() != COMPACT_SHA:
        raise ValueError('Compact arithmetic binding changed')
    spec = importlib.util.spec_from_file_location('pinned_compact_amplification', COMPACT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    tail = module.arithmetic_module()
    return tail.Arithmetic(precision=80), tail


def amplify(tails, direction, kappa_upper, alpha_m, *, number_of_waters=3980,
            retained=None, dual_norms=None, spectral=None):
    """Exact decimal strings are caller-supplied scalar bounds.

    tails: T0,T1,V0,V1. direction: coordinate or strain.
    retained (optional): s,D,W. dual_norms (optional): u_m,v_m,w_water_max,K_abs.
    spectral (optional, requires dual_norms): normalized_lower,alpha_max,alpha_o_max.
    No claim of spectral qualification is made or inferred by this function.
    """
    ar, tail = arithmetic()
    def pos(v, name, strict=False):
        return ar.nonnegative(v, name, positive=strict)
    if type(number_of_waters) is not int or not 1 <= number_of_waters <= 3980:
        raise ValueError('Population outside declared domain')
    if direction not in ('coordinate', 'strain'):
        raise ValueError('Unsupported directional norm')
    if set(tails) != {'T0', 'T1', 'V0', 'V1'}:
        raise ValueError('Exactly four scalar upper bounds required')
    ts = {k: pos(v, k) for k, v in tails.items()}
    T = ts['T0']
    DD, WD = (ts['T1'], T) if direction == 'coordinate' else (ts['V1'], ts['V0'])
    alpha = pos(alpha_m, 'alpha_m', True)
    kappa = pos(kappa_upper, 'kappa_upper', True)
    K0 = ar.rational_power_positive(kappa, 2, 3)
    K1 = ar.mul(ar.integer(5), kappa)
    S = ar.ln(ar.integer(4))
    L = ar.div(ar.mul(ar.integer(12), S), alpha)
    H = ar.div(ar.mul(ar.integer(4), ar.add(ar.mul(ar.integer(14), S),
                                         ar.mul(ar.integer(9), ar.power(S, 2)))), alpha)
    if retained is None:
        sc = ar.mul(ar.integer(number_of_waters), tail.normalization(ar))
        motion = ar.integer(1 if direction == 'coordinate' else 21)
        Dc = ar.mul(ar.mul(ar.div(ar.integer(47), ar.integer(8)), sc), motion)
        Wc = ar.mul(sc, motion)
        retained_source = 'universal_population_contact_weight_upper_bound'
    else:
        if set(retained) != {'s', 'D', 'W'}:
            raise ValueError('Retained upper bounds s,D,W required')
        sc, Dc, Wc = [pos(retained[k], k) for k in ('s', 'D', 'W')]
        retained_source = 'caller_supplied_upper_bounds_not_replayed'
    B0d, B0c = ar.mul(K0, T), ar.mul(K0, sc)
    B1d = ar.add(ar.mul(K0, DD), ar.mul(K1, WD))
    B1c = ar.add(ar.mul(K0, Dc), ar.mul(K1, Wc))
    Eh = ar.mul(ar.integer(3), ar.add(ar.mul(ar.integer(2), B0d), ar.mul(T, B0c)))
    Edh = ar.mul(ar.integer(3), ar.total([
        ar.mul(ar.integer(2), B1d), ar.mul(T, B1c),
        ar.mul(B0d, ar.add(Dc, DD)),
        ar.mul(B0c, ar.add(ar.mul(ar.integer(2), ar.mul(T, Dc)), DD))]))
    Eb = ar.mul(L, T)
    Edb = ar.add(ar.mul(H, ar.mul(T, Dc)), ar.mul(L, DD))
    etaH = ar.mul(alpha, Eb)
    compact_dh = ar.mul(ar.integer(3), ar.add(ar.mul(ar.integer(2), B1c), ar.mul(Dc, B0c)))
    compact_db = ar.mul(L, Dc)
    result = {'status': 'conditional_scalar_amplification', 'direction': direction,
              'universal_parameter_box': True, 'retained_input_source': retained_source,
              'expressions': {k: v.json() for k, v in {
                  'E_h': Eh, 'E_dh': Edh, 'E_beta': Eb, 'E_dbeta': Edb,
                  'E_normalized_H': etaH, 'compact_dh_upper': compact_dh,
                  'compact_dbeta_upper': compact_db, 's_c': sc, 'D_c': Dc,
                  'W_c': Wc, 'V_norm_constant': K0, 'JV_norm_constant': K1,
                  'beta_Lipschitz': L, 'beta_Hessian': H}.items()},
              'explicit_gradient': {'status': 'not_assessed'},
              'induced_dipoles': {'status': 'not_assessed'},
              'total_self_consistent_native_gradient': 'not_assessed'}
    if spectral is not None and dual_norms is None:
        raise ValueError('Spectral amplification requires compact dipole norm bounds')
    if dual_norms is None:
        return result
    if set(dual_norms) != {'u_m', 'v_m', 'w_water_max', 'K_abs'}:
        raise ValueError('Complete explicit-gradient norm input required')
    U, V, W, K = [pos(dual_norms[k], k) for k in ('u_m', 'v_m', 'w_water_max', 'K_abs')]
    fac = ar.div(K, ar.integer(2))
    explicit = ar.mul(fac, ar.add(ar.mul(ar.mul(U, V), Edb), ar.mul(W, Edh)))
    result['explicit_gradient'] = {'status': 'conditional_on_supplied_fixed_dual_norms',
                                   'error_expression': explicit.json()}
    if spectral is None:
        return result
    if set(spectral) != {'normalized_lower', 'alpha_max', 'alpha_o_max'}:
        raise ValueError('Complete separately qualified spectral input required')
    ell, amax, ao = [pos(spectral[k], k, True)
                     for k in ('normalized_lower', 'alpha_max', 'alpha_o_max')]
    if amax.lo < max(alpha.hi, ao.hi):
        raise ValueError('alpha_max does not dominate supplied active values')
    denominator = ar.sub(ell, etaH)
    result['induced_dipoles'] = {'status': 'inconclusive',
                               'qualification_of_lower_bound': 'caller_responsibility_not_replayed',
                               'normalized_lower_after_perturbation': denominator.json()}
    if denominator.lo <= 0:
        return result
    def error(moment):
        numerator = ar.add(ar.mul(ar.sqrt(ao), Eh),
                           ar.mul(ar.sqrt(alpha), ar.mul(Eb, moment)))
        return ar.mul(ar.sqrt(amax), ar.div(numerator, denominator))
    eu, ev = error(U), error(V)
    changed_cross = ar.total([ar.mul(U, ev), ar.mul(V, eu), ar.mul(eu, ev)])
    changed_explicit = ar.mul(fac, ar.add(
        ar.mul(ar.add(compact_db, Edb), changed_cross),
        ar.mul(ar.add(compact_dh, Edh), ar.add(eu, ev))))
    result['induced_dipoles'].update(status='conditional_perturbation_bound',
                                     E_u=eu.json(), E_v=ev.json())
    result['explicit_gradient']['with_changed_duals'] = ar.add(explicit, changed_explicit).json()
    return result
