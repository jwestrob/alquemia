"""Pure fixed Stage A validation arithmetic; no I/O, engines or fitting."""
import math

METALS = ('La', 'Dy')
FAMILIES = ('aquo8', 'hydrated_acetate_mono', 'hydrated_acetate_bi', 'hydrated_acetamide')
CONFIGURATIONS = tuple('validation_' + f + '_orientation_' + side
                       for f in FAMILIES for side in ('minus', 'plus'))
CELLS = tuple(q + '_' + m for q in CONFIGURATIONS for m in METALS)
MISSING = ('unavailable', 'not_attempted', 'failed')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def vector(x, n):
    return isinstance(x, (list, tuple)) and len(x) == n and all(finite(v) for v in x)


def admitted(value, atoms, label):
    if value is None:
        return False
    require(isinstance(value, dict), label + ': record must be a dictionary or None')
    status = value.get('status', 'complete')
    require(status == 'complete' or status in MISSING, label + ': unsupported status')
    if status in MISSING:
        return False
    require(finite(value.get('energy')), label + ': energy absent/nonfinite')
    if atoms is not None:
        force = value.get('forces')
        require(isinstance(force, (list, tuple)) and len(force) == atoms
                and all(vector(v, 3) for v in force), label + ': force shape/nonfinite')
    return True


def check_geometry(data, model, body):
    rows = data['rows']
    require(set(rows) == set(model) == set(CELLS), 'All16 declared data/model slots are required')
    for q in CONFIGURATIONS:
        n = 25 if 'aquo8' in q else 31 if 'acetamide' in q else 26 if 'acetate_bi' in q else 29
        nw = 8 if 'aquo8' in q else 6 if 'acetate_bi' in q else 7
        charge = 2 if 'acetate_' in q else 3
        for metal in METALS:
            gid = q + '_' + metal
            row = rows[gid]
            g = row['geometry']
            require(g['geometry_id'] == gid and g['metal'] == metal
                    and g['atoms'] == n and g['metal_index'] == 0
                    and g['charge'] == charge and g['physical_multiplicity'] == (1 if metal == 'La' else 6)
                    and g['water_count'] == nw, gid + ': chemical state/identity differs')
            require(len(g['positions_A']) == n and all(vector(x, 3) for x in g['positions_A'])
                    and len(g['symbols']) == n and g['symbols'][0] == metal,
                    gid + ': coordinates/symbols invalid')
            bonds = row['input']['bonds_zero_based']
            require(len(bonds) == n and all(isinstance(v, (list, tuple)) for v in bonds), gid + ': graph shape')
            require(all(isinstance(j, int) and not isinstance(j, bool) and 0 <= j < n
                        and i != j and i in bonds[j]
                        for i, ids in enumerate(bonds) for j in ids), gid + ': invalid/asymmetric graph')
            require(all(len(ids) == len(set(ids)) for ids in bonds) and not bonds[0], gid + ': graph duplicates/metal bond')
            require(row['bodies'] == body.bodies(row), gid + ': body mapping differs')
            categories = [b['category'] for b in row['bodies']]
            require(categories.count('ion') == 1 and categories.count('water') == nw
                    and categories.count('ligand') == (0 if 'aquo8' in q else 1), gid + ': body denominator differs')
        a, b = (rows[q + '_' + m] for m in METALS)
        require(a['geometry']['positions_A'] == b['geometry']['positions_A']
                and a['geometry']['symbols'][1:] == b['geometry']['symbols'][1:]
                and a['input']['bonds_zero_based'] == b['input']['bonds_zero_based']
                and a['bodies'] == b['bodies'], q + ': paired geometry/graph/body mismatch')
    for family in FAMILIES:
        for metal in METALS:
            a, b = (rows['validation_' + family + '_orientation_' + s + '_' + metal] for s in ('minus', 'plus'))
            require(a['geometry']['symbols'] == b['geometry']['symbols']
                    and a['input']['bonds_zero_based'] == b['input']['bonds_zero_based']
                    and a['bodies'] == b['bodies'], family + ': orientation atom/body order differs')


def unavailable(dependencies):
    return {'status': 'unavailable', 'passed': None, 'unavailable_dependencies': list(dependencies)}


def energy_record(model_value, reference_value):
    error = model_value - reference_value
    return {'status': 'complete', 'model_kcal_mol': model_value,
            'reference_kcal_mol': reference_value, 'error_kcal_mol': error,
            'passed': abs(error) <= 1.0}


def evaluate(data, model, body, anchors):
    """Exact16 geometry rows and explicit aquo8-central model/reference La/Dy anchors.

    Caller admits provenance and converts E to kcal/mol, F to kcal/mol/Angstrom.
    Reference/model values may be None or explicitly unavailable; no slot is dropped.
    """
    check_geometry(data, model, body)
    require(set(anchors) == {'model', 'reference'}
            and all(set(anchors[k]) == set(METALS) for k in anchors),
            'Explicit model/reference La/Dy aquo8-central anchor slots required')
    rows = data['rows']
    ref = {gid: row.get('reference') for gid, row in rows.items()}
    values = {'model': model, 'reference': ref}
    availability, anchor_availability = {}, {}
    for kind in values:
        for gid in CELLS:
            availability[kind + ':' + gid] = admitted(values[kind][gid], rows[gid]['geometry']['atoms'], kind + ':' + gid)
        for metal in METALS:
            anchor_availability[kind + ':' + metal] = admitted(anchors[kind][metal], None, kind + ' aquo8:' + metal)
    cells = {gid: {'status': 'complete' if all(availability[k + ':' + gid] for k in values) else 'unavailable',
                   'model_available': availability['model:' + gid], 'reference_available': availability['reference:' + gid]}
             for gid in CELLS}

    def missing(gids, use_anchors=False):
        out = [kind + ':' + gid for kind in values for gid in gids if not availability[kind + ':' + gid]]
        if use_anchors:
            out += ['aquo8_central:' + k for k, ok in anchor_availability.items() if not ok]
        return out

    def delta(kind, q):
        return values[kind][q + '_Dy']['energy'] - values[kind][q + '_La']['energy']

    exchange = {}
    for q in CONFIGURATIONS:
        absent = missing([q + '_' + m for m in METALS], True)
        if absent:
            exchange[q] = unavailable(absent)
        else:
            c = {k: delta(k, q) - (anchors[k]['Dy']['energy'] - anchors[k]['La']['energy']) for k in values}
            exchange[q] = energy_record(c['model'], c['reference'])

    work = {}
    for family in FAMILIES:
        prefix = 'validation_' + family + '_orientation_'
        for kind in ('La', 'Dy', 'paired'):
            metals = METALS if kind == 'paired' else (kind,)
            gids = [prefix + s + '_' + m for m in metals for s in ('minus', 'plus')]
            absent = missing(gids)
            key = family + '_' + kind
            if absent:
                work[key] = unavailable(absent)
            else:
                def w(which, metal):
                    return values[which][prefix + 'plus_' + metal]['energy'] - values[which][prefix + 'minus_' + metal]['energy']
                v = {which: (w(which, 'Dy') - w(which, 'La') if kind == 'paired' else w(which, kind)) for which in values}
                work[key] = energy_record(v['model'], v['reference'])

    per = {k: {} for k in ('La', 'Dy', 'paired')}
    atomic = {}
    for kind in per:
        for q in CONFIGURATIONS:
            gids = [q + '_' + m for m in (METALS if kind == 'paired' else (kind,))]
            absent = missing(gids)
            if absent:
                per[kind][q] = unavailable(absent)
                atomic[q + '_' + kind] = unavailable(absent)
                continue
            row = rows[q + '_La']  # Matched positions/body mapping checked for both metals.
            if kind == 'paired':
                mf = body.subtract(model[q + '_Dy']['forces'], model[q + '_La']['forces'])
                rf = body.subtract(ref[q + '_Dy']['forces'], ref[q + '_La']['forces'])
            else:
                mf, rf = model[q + '_' + kind]['forces'], ref[q + '_' + kind]['forces']
            weighted, raw = body.rigid_residual(row, mf, rf)
            rms = 2 * math.sqrt(sum(v * v for v in weighted))
            maximum = max(abs(v) for v in raw)
            per[kind][q] = {'status': 'complete', 'rigid_RMS_kcal_mol_A': rms,
                'maximum_body_component_error_kcal_mol_A': maximum, 'passed': rms <= 2 and maximum <= 5,
                'raw_body_components_kcal_mol_A': raw, 'model_body_response': body.response(row, mf),
                'reference_body_response': body.response(row, rf)}
            residual = body.subtract(mf, rf)
            atomic[q + '_' + kind] = {'status': 'complete', 'atomic_force_residual_kcal_mol_A': residual,
                'atom_vector_RMS_kcal_mol_A': math.sqrt(sum(x*x for v in residual for x in v) / len(residual)),
                'used_as_pass_criterion': False}

    groups = {'balanced_exchange': list(exchange.values()), 'orientation_work': list(work.values()),
              'rigid_response': [v for group in per.values() for v in group.values()]}
    counts = {name: {'declared': len(group), 'available': sum(v['status'] == 'complete' for v in group),
                     'passed': sum(v['passed'] is True for v in group), 'failed': sum(v['passed'] is False for v in group)}
              for name, group in groups.items()}
    complete = all(v['status'] == 'complete' for v in cells.values()) and all(anchor_availability.values())
    passed = complete and all(v['passed'] is True for group in groups.values() for v in group)
    return {'status': 'complete' if complete else 'partial', 'validation_passed': passed,
            'calculated_failure_observed': any(v['passed'] is False for group in groups.values() for v in group),
            'cells': cells, 'anchor_availability': anchor_availability, 'balanced_exchange': exchange,
            'orientation_work': work, 'per_fixture': per, 'atomic_diagnostics': atomic, 'counts': counts,
            'denominators': {'cells': 16, 'configurations': 8, 'balanced_exchange': 8, 'orientation_work': 12,
                             'rigid_response': 24, 'total_comparisons': 44},
            'limits': {'energy_absolute_kcal_mol': 1.0, 'rigid_RMS_kcal_mol_A': 2.0,
                       'maximum_body_component_kcal_mol_A': 5.0},
            'reference_scope': 'Conditional on admitted frozen-f ECP-PBE0-D4 energies/forces and their retained method/state limits',
            'fit_objective_defined': False, 'physical_accuracy_qualified': False, 'stage_B_released': False,
            'production_released': False}
