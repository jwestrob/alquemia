"""Deterministic complete-poll search. Native work is supplied by a bounded callback."""
import concurrent.futures
import hashlib
import json
import math

VECTOR_CAP, SEARCH_CAP = 241, 233


def token(values):
    return hashlib.sha256(json.dumps(values, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def eigenvalues(matrix):
    """Small symmetric Jacobi eigensolver for the normalized sensitivity Gram matrix."""
    scale = max(abs(v) for r in matrix for v in r)
    if scale == 0:
        return [0.]*len(matrix)
    a = [[v/scale for v in r] for r in matrix]
    n = len(a)
    for _ in range(100*n*n):
        i, j = max(((i, j) for i in range(n) for j in range(i+1, n)), key=lambda z: abs(a[z[0]][z[1]]))
        if abs(a[i][j]) <= 1e-15:
            break
        angle = .5*math.atan2(2*a[i][j], a[j][j]-a[i][i])
        c, s = math.cos(angle), math.sin(angle)
        x, y, z = a[i][i], a[j][j], a[i][j]
        a[i][i] = c*c*x-2*c*s*z+s*s*y
        a[j][j] = s*s*x+2*c*s*z+c*c*y
        a[i][j] = a[j][i] = 0.
        for k in range(n):
            if k in (i, j):
                continue
            x, y = a[k][i], a[k][j]
            a[k][i] = a[i][k] = c*x-s*y
            a[k][j] = a[j][k] = s*x+c*y
    else:
        raise ValueError('Sensitivity eigensolver did not converge')
    return sorted((max(0., a[i][i]*scale) for i in range(n)), reverse=True)


def search(parameters, evaluate, stop_event, update, workers=8):
    """evaluate(values) returns {'score': ..., 'model': ...}; all callbacks finish before selection."""
    if len(parameters) != 4 or not 1 <= workers <= 8:
        raise ValueError('Four parameters and at most eight workers required')
    initial = [p['initial'] for p in parameters]
    z0 = [(p['initial']-p['lower'])/(p['upper']-p['lower']) for p in parameters]
    if not all(math.isfinite(v) and 0 <= v <= 1 for v in z0):
        raise ValueError('Initial parameter outside bounds')
    records = {}
    state = {'status': 'searching', 'vector_cap': VECTOR_CAP, 'search_cap': SEARCH_CAP,
             'evaluations': records, 'polls': [], 'candidate_frozen': False,
             'numerical_optimization_complete': False, 'physical_accuracy_qualified': False}

    def values(z):
        # Preserve exact baseline values rather than round-tripping interval arithmetic.
        if z == z0:
            return initial[:]
        return [p['lower']+v*(p['upper']-p['lower']) for p, v in zip(parameters, z)]

    def identity(z):
        return token(values(z))

    def batch(points, role):
        todo, ordered = [], []
        for z in points:
            key = identity(z)
            ordered.append(key)
            if key in records:
                if records[key]['status'] != 'complete':
                    raise ValueError('Attempt to repeat unavailable vector')
                continue
            records[key] = {'parameters': values(z), 'normalized_parameters': z[:], 'role': role,
                            'status': 'scheduled', 'model': None, 'score': None}
            todo.append(key)
        if len(records) > VECTOR_CAP:
            raise ValueError('Total vector cap exceeded')
        update(state)
        errors = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            def call(key):
                if stop_event.is_set():
                    raise RuntimeError('Native scheduling stopped before vector')
                return evaluate(records[key]['parameters'])
            futures = {pool.submit(call, key): key for key in todo}
            for future in concurrent.futures.as_completed(futures):
                key = futures[future]
                try:
                    value = future.result()
                    if not math.isfinite(value['score']['objective']):
                        raise ValueError('Nonfinite score')
                    records[key].update(status='complete', **value)
                except Exception as exc:
                    stop_event.set()
                    records[key].update(status='unavailable', error=repr(exc))
                    errors.append((key, repr(exc)))
                    for other in futures:
                        if other is not future:
                            other.cancel()
                update(state)
        if errors:
            state.update(status='failed_or_unavailable', errors=errors)
            update(state)
            raise RuntimeError('Vector/native failure; no further polls')
        return [records[key] for key in ordered]

    best = batch([z0], 'baseline')[0]
    state['initial_key'] = identity(z0)
    step, stop = .10, 'budget_limited'
    while True:
        center = best['normalized_parameters']
        points, seen = [], set()
        for j in range(4):
            for sign in (-1, 1):
                z = center[:]
                z[j] = min(1., max(0., z[j]+sign*step))
                key = identity(z)
                if key == token(best['parameters']) or key in seen:
                    continue
                seen.add(key)
                points.append(z)
        required = sum(identity(z) not in records for z in points)
        if len(records)+required > SEARCH_CAP:
            state['unstarted_poll'] = {'step': step, 'points': points, 'new_vectors_required': required,
                                       'reason': 'complete_poll_exceeds_remaining_search_budget'}
            break
        if not points:
            stop = 'no_distinct_poll_points'
            break
        origin = token(best['parameters'])
        trials = batch(points, 'search')
        chosen = best
        for rec in trials:  # Original poll order, independent of completion order.
            if rec['score']['objective'] < chosen['score']['objective']:
                chosen = rec
        improved = chosen is not best
        best = chosen
        state['polls'].append({'step': step, 'complete': True, 'origin': origin,
            'trial_keys': [token(r['parameters']) for r in trials], 'selected': token(best['parameters']),
            'improved': improved, 'best_objective': best['score']['objective']})
        state['best_key'] = token(best['parameters'])
        update(state)
        if not improved:
            step /= 2
        if step < .00625:
            stop = 'local_step_rule_met'
            break
    selected = token(best['parameters'])
    state.update(status='sensitivity', search_stopping_rule=stop, best_key=selected,
                 search_distinct_vectors=len(records), numerical_optimization_complete=stop == 'local_step_rule_met')
    perturbations = []
    for j in range(4):
        lo, hi = best['normalized_parameters'][:], best['normalized_parameters'][:]
        lo[j], hi[j] = max(0., lo[j]-.02), min(1., hi[j]+.02)
        perturbations.extend((lo, hi))
    sampled = batch(perturbations, 'sensitivity')
    columns = []
    for j in range(4):
        a, b = sampled[2*j]['score']['weighted_residuals'], sampled[2*j+1]['score']['weighted_residuals']
        width = perturbations[2*j+1][j]-perturbations[2*j][j]
        if width <= 0 or len(a) != len(b):
            raise ValueError('Invalid sensitivity displacement')
        columns.append([(y-x)/width for x, y in zip(a, b)])
    gram = [[sum(a*b for a, b in zip(x, y)) for y in columns] for x in columns]
    singular = [math.sqrt(v) for v in eigenvalues(gram)]
    rank = sum(v > singular[0]*1e-6 for v in singular) if singular[0] else 0
    condition = singular[0]/singular[-1] if singular[-1] else None
    near_bound = [min(z, 1-z) <= .01 for z in best['normalized_parameters']]
    state.update(status='finite_pilot_complete', best_key=selected, distinct_vectors=len(records),
        candidate_parameters=best['parameters'], candidate_frozen=True, near_bound_by_parameter=near_bound,
        optimizer_budget_limited=stop == 'budget_limited',
        sensitivity={'singular_values': singular, 'rank': rank, 'condition': condition,
            'column_norms': [math.sqrt(sum(v*v for v in c)) for c in columns],
            'weakly_identified': rank < 4 or condition is None or condition > 1e6,
            'parameter_intervals_used': [[perturbations[2*j][j], perturbations[2*j+1][j]] for j in range(4)],
            'candidate_not_reselected_from_sensitivity': True},
        development_checks_passed=best['score']['development_checks_passed'],
        identifiability_and_bounds_passed=rank == 4 and condition is not None and condition <= 1e6 and not any(near_bound),
        automatic_validation_release=False, production_released=False)
    update(state)
    return state
