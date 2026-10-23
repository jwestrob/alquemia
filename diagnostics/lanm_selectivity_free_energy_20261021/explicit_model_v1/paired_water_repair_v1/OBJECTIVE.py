"""Four-cell paired water objective. Saved arithmetic only; no native engines."""
import math

GEOMETRIES = ('water_central', 'aquo8_central')
CELLS = tuple(q+'_'+m for q in GEOMETRIES for m in ('La', 'Dy'))
WEIGHTS = {'exchange': .25, 'La': .25, 'Dy': .25, 'paired': .25}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def admit(data):
    rows = data['rows']
    for q in GEOMETRIES:
        a, b = rows[q+'_La'], rows[q+'_Dy']
        require(a['geometry']['positions_A'] == b['geometry']['positions_A'], 'Unmatched coordinates')
        require(a['input']['bonds_zero_based'] == b['input']['bonds_zero_based'], 'Unmatched bonds')
        require(a['bodies'] == b['bodies'], 'Unmatched rigid bodies')
        for row, metal in ((a, 'La'), (b, 'Dy')):
            g = row['geometry']
            require(g['metal'] == metal and g['charge'] == 3, 'Water state changed')
            require(g['atoms'] == (4 if q == 'water_central' else 25), 'Water atom count')
            require(g['symbols'][g['metal_index']] == metal, 'Metal index')
            require(all(bd['category'] in ('ion', 'water') for bd in row['bodies']), 'Non-water body')
            require(sorted(i for bd in row['bodies'] for i in bd['indices']) == list(range(g['atoms'])),
                    'Incomplete body coverage')
            for source in ('reference', 'baseline'):
                require(math.isfinite(row[source]['energy']), 'Nonfinite energy')
                check_forces(row, row[source]['forces'])
    return rows


def check_forces(row, forces):
    require(len(forces) == row['geometry']['atoms'] and
            all(len(v) == 3 and all(math.isfinite(x) for x in v) for v in forces), 'Incomplete forces')


def delta(values, q):
    return values[q+'_Dy']['energy'] - values[q+'_La']['energy']


def contrast(values):
    return delta(values, 'water_central') - delta(values, 'aquo8_central')


def evaluate(data, model, body_helper):
    rows = admit(data)
    require(set(model) == set(CELLS), 'Exactly four model cells required')
    reference = {k: rows[k]['reference'] for k in CELLS}
    for k in CELLS:
        require(math.isfinite(model[k]['energy']), 'Nonfinite model energy')
        check_forces(rows[k], model[k]['forces'])
    ce = contrast(model) - contrast(reference)
    vectors = {'exchange': [ce]}
    controls, per_fixture, internal = {}, {}, {}
    for kind in ('La', 'Dy', 'paired'):
        vec, raw = [], []
        per_fixture[kind] = {}
        for q in GEOMETRIES:
            row = rows[q+'_La']
            if kind == 'paired':
                mf = body_helper.subtract(model[q+'_Dy']['forces'], model[q+'_La']['forces'])
                rf = body_helper.subtract(reference[q+'_Dy']['forces'], reference[q+'_La']['forces'])
            else:
                mf, rf = model[q+'_'+kind]['forces'], reference[q+'_'+kind]['forces']
            weighted, components = body_helper.rigid_residual(row, mf, rf)
            vec.extend(v / math.sqrt(len(GEOMETRIES)) for v in weighted)
            raw.extend(components)
            per_fixture[kind][q] = {'rigid_RMS_kcal_mol_A': 2*math.sqrt(sum(v*v for v in weighted)),
                'maximum_body_component_error_kcal_mol_A': max(abs(v) for v in components),
                'raw_body_components_kcal_mol_A': components,
                'model_body_response': body_helper.response(row, mf),
                'reference_body_response': body_helper.response(row, rf)}
            residual = body_helper.subtract(mf, rf)
            internal[q+'_'+kind] = {'atomic_force_residual_kcal_mol_A': residual,
                'atom_vector_RMS_kcal_mol_A': math.sqrt(sum(x*x for v in residual for x in v)/len(residual)),
                'used_as_internal_force_target': False}
        vectors[kind] = vec
        rms, maximum = 2*math.sqrt(sum(v*v for v in vec)), max(abs(v) for v in raw)
        controls[kind] = {'category_balanced_rigid_RMS_kcal_mol_A': rms,
            'max_body_component_error_kcal_mol_A': maximum, 'passed': rms <= 2 and maximum <= 5}
    weighted = [math.sqrt(WEIGHTS[k])*v for k in WEIGHTS for v in vectors[k]]
    require(all(math.isfinite(v) for v in weighted), 'Nonfinite objective')
    return {'objective': sum(v*v for v in weighted), 'weighted_residuals': weighted,
        'blocks': {k: sum(v*v for v in vectors[k]) for k in WEIGHTS}, 'block_weights': WEIGHTS,
        'balanced_exchange_error_kcal_mol': ce, 'model_balanced_exchange_kcal_mol': contrast(model),
        'reference_balanced_exchange_kcal_mol': contrast(reference), 'controls': controls,
        'per_fixture': per_fixture, 'atomic_diagnostics': internal,
        'development_checks_passed': abs(ce) <= 1 and all(c['passed'] for c in controls.values()),
        'physical_accuracy_qualified': False, 'production_released': False}


def method_sensitivity(data, model, water48, panel):
    """Use exact Water48 outward axis; do not use the inward training axis."""
    rows = admit(data)
    projected, baseline, refs = {}, {}, {}
    for metal in ('La', 'Dy'):
        row = rows['water_central_'+metal]
        central = panel['central'][metal]
        require(central['origin_cell']['positions_A'] == row['geometry']['positions_A'], 'Water48 center geometry')
        require(central['origin_cell']['symbols'] == row['geometry']['symbols'], 'Water48 atom identity')
        u = central['unit_metal_to_O']
        water = [b for b in row['bodies'] if b['category'] == 'water']
        require(len(water) == 1, 'Dimer water body')
        ids = water[0]['indices']
        def slope(forces):
            return -sum(sum(forces[i][k] for i in ids)*u[k] for k in range(3))
        projected[metal] = slope(model['water_central_'+metal]['forces'])
        baseline[metal] = slope(row['baseline']['forces'])
        saved = water48['per_metal'][metal]
        require(abs(baseline[metal]-saved['central_AMOEBA_analytic']['coordinate_derivative_kcal_mol_A']) <= 1e-10,
                'Saved AMOEBA directional replay')
        require(saved['methods']['pbe0']['sampled_analytic_consistency_at_both_steps'] and
                all(all(saved['methods'][m]['two_step_resolved'].values()) for m in ('pbe0', 'cc')), 'Water48 consistency')
        refs[metal] = {label: saved['methods'][method]['prespecified_smaller_step_summary'][key]
            for label, method, key in [('PBE0_D4_uncorrected', 'pbe0', 'complex_D4'),
                ('PBE0_D4_CP', 'pbe0', 'CP_D4'), ('CC_uncorrected', 'cc', 'complex'), ('CC_CP', 'cc', 'CP')]}
    projected['paired'] = projected['Dy']-projected['La']
    baseline['paired'] = baseline['Dy']-baseline['La']
    refs['paired'] = {k: refs['Dy'][k]-refs['La'][k] for k in refs['La']}
    result = {}
    for kind in ('La', 'Dy', 'paired'):
        result[kind] = {'model_dE_dq_kcal_mol_A': projected[kind], 'baseline_dE_dq_kcal_mol_A': baseline[kind],
            'reference_dE_dq_kcal_mol_A': refs[kind], 'comparisons': {}}
        for method, target in refs[kind].items():
            residual, original = projected[kind]-target, baseline[kind]-target
            result[kind]['comparisons'][method] = {'residual_kcal_mol_A': residual,
                'baseline_residual_kcal_mol_A': original, 'wrong_sign': projected[kind]*target < 0,
                'absolute_residual_worsened': abs(residual) > abs(original)}
    return {'outward_water_translation_positive': True, 'results': result,
        'CC_residual_deterioration_flag': any(v['absolute_residual_worsened'] for r in result.values()
            for k, v in r['comparisons'].items() if k.startswith('CC')),
        'used_in_objective_or_candidate_selection': False, 'method_accuracy_qualified': False}
