"""Complete, fixed-denominator pair/cooperative contrasts in the input unit."""
import math

GIDS = tuple(f'aquo{n}_central_{m}' for n in (8, 9) for m in ('La', 'Dy'))
SIGNS = dict(zip(GIDS, (1, -1, -1, 1)))


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def decompose(clusters, pairs, ions):
    if set(clusters) != set(GIDS) or set(pairs) != set(GIDS) or set(ions) != {'La', 'Dy'}:
        raise ValueError('Incomplete or extra source/ion inventory')
    for gid in GIDS:
        n = int(gid[4])
        if set(pairs[gid]) != set(range(n)):
            raise ValueError('Incomplete or extra water inventory: ' + gid)
    values = list(clusters.values()) + list(ions.values())
    values += [v for gid in GIDS for v in pairs[gid].values()]
    if not all(finite(x) for x in values):
        raise ValueError('Nonfinite or invalid energy')
    sums = {g: math.fsum(pairs[g][i] for i in range(int(g[4]))) for g in GIDS}
    pair_delta = math.fsum(SIGNS[g] * sums[g] for g in GIDS)
    ion_delta = ions['Dy'] - ions['La']
    c = math.fsum(SIGNS[g] * clusters[g] for g in GIDS)
    p = pair_delta - ion_delta
    b = c - p
    return {'cluster_contrast': c, 'pair_contrast': p, 'cooperative_remainder': b,
            'pair_energy_sums': sums, 'pair_sum_contrast_before_ion_subtraction': pair_delta,
            'own_ion_difference_subtracted': ion_delta,
            'closure_residual': c - math.fsum((p, b))}


def compare(reference, model):
    keys = ('cluster_contrast', 'pair_contrast', 'cooperative_remainder')
    if not all(finite(x[k]) for x in (reference, model) for k in keys):
        raise ValueError('Unavailable decomposition')
    errors = {k: model[k] - reference[k] for k in keys}
    errors['closure_residual'] = errors['cluster_contrast'] - math.fsum(
        (errors['pair_contrast'], errors['cooperative_remainder']))
    return errors
