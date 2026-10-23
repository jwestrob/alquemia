"""Declared joint engineering objective; pure saved-value arithmetic, no engines."""
import copy
import math

PARENTS = ('aquo8_central_La', 'aquo8_central_Dy',
           'aquo9_central_La', 'aquo9_central_Dy')
STATES = ('La_dimer', 'La_trimer', 'Dy_dimer', 'Dy_trimer')
BLOCK_WEIGHTS = {'old_EF': .25, 'pair': .25, 'C': .25, 'chi': .25}
RESIDUAL_COUNT = 531


def require(ok, message):
    if not ok:
        raise ValueError(message)


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def vector(x):
    return isinstance(x, (list, tuple)) and len(x) == 3 and all(number(v) for v in x)


def targets_from_inventory(inventory):
    """Extract the already-pinned target inventory; no estimator selection by outcome."""
    groups = inventory['groups']
    source = groups['pair_and_collective_dipoles']
    require(set(source['parents']) == set(PARENTS), 'Exactly four parent contexts')
    pairs = [{'id': r['id'], 'parent': r['parent'], 'water_index': r['water_index'],
              'reference_Debye': copy.deepcopy(r['reference_pair_mu_Debye'])}
             for r in source['pair_cycles']]
    require(len(pairs) == 34 and len({r['id'] for r in pairs}) == 34, 'Exactly34 literal pairs')
    for parent in PARENTS:
        require(sum(r['parent'] == parent for r in pairs) == (8 if parent.startswith('aquo8_') else 9),
                'Complete literal pairs within each parent')
    require(all(r['parent'] in PARENTS and vector(r['reference_Debye']) for r in pairs),
            'Finite pair reference vectors')
    C = {p: copy.deepcopy(source['parents'][p]['C']['reference_mu_Debye']) for p in PARENTS}
    require(all(vector(v) for v in C.values()), 'Finite full-shell C references')
    fields = copy.deepcopy(groups['conditional_susceptibility']['reported_rows'])
    quantities = set(STATES) | {'La', 'Dy', 'Dy_minus_La'}
    require(len(fields) == 28 and
            {(r['quantity'], r['step'], r['estimator']) for r in fields} ==
            {(q, s, e) for q in quantities for s in ('h', 'half') for e in ('energy', 'moment')},
            'Exactly28 original field rows')
    require(all(number(r['reference_a0_cubed']) and
                r['field_size_au'] == (.001 if r['step'] == 'h' else .0005) for r in fields),
            'Finite field references and original steps')
    selected = [r for r in fields if r['quantity'] in STATES and r['step'] == 'h'
                and r['field_size_au'] == .001 and r['estimator'] == 'energy' and r['kind'] == 'state']
    require(len(selected) == 4 and {r['quantity'] for r in selected} == set(STATES),
            'Four fixed h=.001 energy-curvature targets')
    return {'pair_rows': pairs, 'C_Debye': C,
            'chi_a0_cubed': {r['quantity']: r['reference_a0_cubed'] for r in selected},
            'field_rows': fields}


def response_blocks(response, targets):
    """Response fields are physical quantities, not preweighted residuals."""
    require(set(response) == {'pair_dipoles_Debye', 'C_dipoles_Debye', 'chi_states_a0_cubed'},
            'Exact response input fields')
    pair = response['pair_dipoles_Debye']; C = response['C_dipoles_Debye']; chi = response['chi_states_a0_cubed']
    require(set(pair) == {r['id'] for r in targets['pair_rows']} and len(pair) == 34,
            'All34 literal pair predictions required')
    require(set(C) == set(PARENTS) and set(chi) == set(STATES), 'All four C and four chi predictions required')
    require(all(vector(v) for v in pair.values()) and all(vector(v) for v in C.values())
            and all(number(v) for v in chi.values()), 'Missing or nonfinite response prediction')
    residuals = {'pair': [], 'C': [], 'chi': []}; details = {'pair': {}, 'C': {}, 'chi': {}}
    contexts = {}
    for parent in PARENTS:
        rows = [r for r in targets['pair_rows'] if r['parent'] == parent]
        errors = []
        for row in rows:
            error = [x-y for x, y in zip(pair[row['id']], row['reference_Debye'])]
            errors.extend(error)
            details['pair'][row['id']] = {'parent': parent, 'model_Debye': list(pair[row['id']]),
                                          'reference_Debye': row['reference_Debye'], 'error_Debye': error}
            residuals['pair'].extend(x/math.sqrt(4*len(rows)*3) for x in error)
        contexts[parent] = sum(x*x for x in errors)/(len(rows)*3)
        error = [x-y for x, y in zip(C[parent], targets['C_Debye'][parent])]
        details['C'][parent] = {'model_Debye': list(C[parent]), 'reference_Debye': targets['C_Debye'][parent],
                                'error_Debye': error}
        residuals['C'].extend(x/math.sqrt(12) for x in error)
    for state in STATES:
        error = chi[state]-targets['chi_a0_cubed'][state]
        details['chi'][state] = {'model_a0_cubed': chi[state], 'reference_a0_cubed': targets['chi_a0_cubed'][state],
                                 'error_a0_cubed': error, 'reference_estimator': 'energy', 'field_size_au': .001}
        residuals['chi'].append(error/2)
    derived = dict(chi)
    for metal in ('La', 'Dy'):
        derived[metal] = chi[metal+'_trimer']-chi[metal+'_dimer']
    derived['Dy_minus_La'] = derived['Dy']-derived['La']
    fields = []
    for row in targets['field_rows']:
        used = row['quantity'] in STATES and row['step'] == 'h' and row['estimator'] == 'energy'
        fields.append(dict(row, candidate_a0_cubed=derived[row['quantity']],
                           candidate_minus_reference_a0_cubed=derived[row['quantity']]-row['reference_a0_cubed'],
                           used_in_objective=used))
    blocks = {k: sum(x*x for x in values) for k, values in residuals.items()}
    require(all(number(v) for v in blocks.values()), 'Finite response block objectives required')
    return {'objectives': blocks, 'residuals': residuals, 'per_observable': details,
            'pair_context_mean_squared_errors': contexts, 'field_rows': fields}


def evaluate(data, model, response, targets, old_objective, body):
    """Delegate the frozen413-vector EF scorer; append only declared response blocks."""
    old = old_objective.evaluate(data, model, body)
    require(len(old['weighted_residuals']) == 413 and old['block_weights'] ==
            {'energy': .25, 'La': .25, 'Dy': .25, 'paired': .25}
            and old['available_energy_controls'] == 10 and old['available_force_controls'] == 15,
            'Unchanged original413-vector and25-gate scorer required')
    r = response_blocks(response, targets)
    blocks = {'old_EF': old['objective'], **r['objectives']}
    require(all(number(v) for v in blocks.values()) and number(sum(blocks.values())),
            'Finite joint block objectives required')
    weighted = [x/2 for x in old['weighted_residuals']]
    for k in ('pair', 'C', 'chi'):
        weighted.extend(x/2 for x in r['residuals'][k])
    require(len(weighted) == RESIDUAL_COUNT and all(number(x) for x in weighted), 'Finite531-vector required')
    return {'objective': sum(blocks.values())/4, 'block_objectives': blocks,
            'block_weights': dict(BLOCK_WEIGHTS), 'weighted_residuals': weighted,
            'old_EF': old, 'response': r,
            'residual_layout': {'old_EF': [0, 413], 'pair': [413, 515], 'C': [515, 527], 'chi': [527, 531],
                                'parent_order': list(PARENTS), 'chi_order': list(STATES),
                                'pair_order': [r['id'] for p in PARENTS for r in targets['pair_rows'] if r['parent'] == p]},
            'engineering_scales': {'pair_Debye': 1, 'C_Debye': 1, 'chi_a0_cubed': 1},
            'original_EF_checks_passed': old['development_checks_passed'],
            'response_physical_acceptance': None, 'joint_physical_acceptance': None,
            'excluded_from_objective': ['P', 'B', 'derived susceptibility increments and paired contrast',
                                         'other24 field rows', 'eight E3/mu3 cycles'],
            'interpretation': 'Engineering weights/scales, not variances, confidence limits or independent-observation likelihood.',
            'production_released': False}
