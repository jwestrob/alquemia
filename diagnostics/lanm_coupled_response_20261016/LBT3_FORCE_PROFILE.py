#!/usr/bin/env python3
"""Profile saved whole-LBT3 forces and prepare one paired torsion derivative check.

Slurm only; no energy engine, optimizer or force-selected displacement proposal.
CLI: python LBT3_FORCE_PROFILE.py --repo-root REPO --config CONFIG --output FRESH
Config pins: plan, modes_result, solvent_result, helpers (the three mode helpers
plus scaffold_accommodation.py). Optional implementation pins this adapter.
Selection is fixed below and must be declared in the pinned plan before running.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import traceback

SOURCES = ('7CCO', '7CCN')
METALS = ('La', 'Dy')
EV_TO_KCAL = 23.06054783061903
RCOND = 1e-10
RESIDUAL_MIN = 1e-5
STEPS_A = (0.001, 0.0005)
RULE = {
    'load': 'dot(force, J_k / ||J_k||_all256); positive force follows positive named coordinate',
    'torsion_score': 'maximum over both sources of abs(La), abs(Dy), abs(Dy-minus-La) composite normalized force loads',
    'tie_break': 'saved mode order',
    'mandatory': ['metal_x', 'metal_y', 'metal_z'],
    'maximum_torsions': 4,
    'independence': 'greedy normalized-column residual norm > 1e-5 in BOTH source metrics',
    'FD_nominee': 'largest absolute composite normalized Dy-minus-La NONMETAL load in 7CCO; saved mode order breaks ties',
    'purpose': 'common candidate coordinate basis; no differential optimization or finite-energy prediction',
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(record):
    path = Path(record['path']).resolve(strict=True)
    require(path.is_file() and pin(path)['sha256'] == record['sha256'], 'Missing/changed input: ' + record['path'])
    return path


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def unique(rows, key):
    result = {key(row): row for row in rows}
    require(len(result) == len(rows), 'Duplicate input rows')
    return result


def metric(matrix, np):
    gram = matrix.T @ matrix
    eigenvalues = np.linalg.eigvalsh(gram)
    keep = eigenvalues > RCOND * eigenvalues[-1]
    require(eigenvalues[0] >= -1e-10 * eigenvalues[-1] and keep.any(), 'Invalid displacement metric')
    return {'normalized_Gram': gram.tolist(), 'Gram_eigenvalues': eigenvalues.tolist(),
            'retained_span_rank': int(keep.sum()), 'relative_cutoff': RCOND,
            'retained_Gram_condition': float(eigenvalues[keep].max() / eigenvalues[keep].min()),
            'interpretation': 'equal-weight all256 Cartesian displacement overlap; not stiffness or energy curvature'}


def source_profile(source, mode_row, prep_row, corrections, force_by, raw_rows, gradient_units, np, Kinematics):
    modes = read(mode_row['modes']); prepared = read(prep_row['preparation'])
    require(mode_row['status'] == 'prepared' and modes['source_id'] == prepared['source_id'] == source
            and modes['source_preparation'] == prep_row['preparation'], 'Mode/source preparation binding differs')
    require(prepared['status'] == 'prepared' and prepared['atoms'] == 256 and prepared['charge'] == -2
            and prepared['protein_atoms'] == 255 and prepared['protein_charge'] == -5
            and prepared['waters'] == prepared['spectators'] == prepared['caps'] == [], 'Whole source chemistry differs')
    atoms = read(prepared['physical_atoms']); native = read(prepared['native_peptide_atoms'])
    bonds = read(prepared['covalent_bonds']); verify(prepared['native_peptide_system'])
    data = read(modes['kinematics']); saved_jac = np.load(verify(modes['origin_jacobian']), allow_pickle=False)
    ids = [a['id'] for a in atoms]; x = np.asarray([a['xyz_A'] for a in atoms])
    require(len(ids) == len(set(ids)) == 256 and ids[-1] == 'metal'
            and ids == modes['physical_source_ids'] == data['physical_ids']
            and [a['id'] for a in native] == ids[:255], 'Physical atom ordering differs')
    require(modes['physical_source'] == prepared['physical_atoms'] and modes['covalent_bonds'] == prepared['covalent_bonds']
            and data['positions_A'] == x.tolist() and data['metal_index'] == 255, 'Physical origin/map differs')
    require(data['modes'] == modes['modes'] and [m['name'] for m in data['modes']] == modes['mode_order']
            and len(data['modes']) == 64 and modes['mode_counts'] ==
            {'phi': 16, 'psi': 16, 'sidechain': 29, 'metal_translation': 3}, 'Catalogue differs')
    require(data['bonds'] == [b['native_particle_indices'] for b in bonds]
            and all([ids[i] for i in b['native_particle_indices']] == [b['atom_a_id'], b['atom_b_id']] for b in bonds),
            'Native bond mapping differs')
    kin = Kinematics(data); replay, _, jac, _ = kin.evaluate(np.zeros(64))
    require(jac.shape == (64, 256, 3) and np.isfinite(jac).all()
            and np.allclose(replay, x, rtol=0, atol=1e-12)
            and np.allclose(jac, saved_jac, rtol=0, atol=1e-12), 'Source/Jacobian replay differs')
    matrix = jac.reshape(64, -1).T; scales = np.linalg.norm(matrix, axis=0); normalized = matrix / scales
    require(np.all(scales > 1e-12) and np.allclose(scales, modes['column_norms'], rtol=0, atol=1e-10)
            and np.allclose(normalized.T @ normalized, modes['normalized_Gram'], rtol=0, atol=1e-10), 'Saved metric differs')
    forces = {component: {} for component in ('composite', 'MACE_vacuum', 'solvent')}
    for metal in METALS:
        task = prepared['endpoints']['whole'][metal]; item = corrections[(source, metal)]
        raw = force_by[task['task_id']]
        require(item['status'] == 'available' and item['MACE_forces_available'] and raw['status'] == 'complete'
                and item['origin_task_id'] == task['task_id'] and item['xyz'] == task['xyz'], 'Unavailable/different force origin')
        require(all(raw[key] == value for key, value in task.items()) and task['charge'] == -2
                and task['metal_index'] == 255 and task['multiplicity'] == (1 if metal == 'La' else 6), 'Force physical state differs')
        lines = verify(task['xyz']).read_text().splitlines(); rows = [line.split() for line in lines[2:] if line.strip()]
        require(int(lines[0]) == len(rows) == 256 and all(len(r) == 4 for r in rows)
                and [r[0] for r in rows] == [metal if a['element'] == 'M' else a['element'] for a in atoms]
                and [[float(v) for v in r[1:]] for r in rows] == x.tolist(), 'Force XYZ/source ordering or coordinates differ')
        mace = np.asarray(item['MACE_forces_eV_A']) * EV_TO_KCAL
        solvent = -np.asarray(item['solvent_gradient_kcal_mol_A'])
        composite = np.asarray(item['composite_forces_eV_A']) * EV_TO_KCAL
        require(all(v.shape == (256, 3) and np.isfinite(v).all() for v in (mace, solvent, composite)), 'Invalid force array')
        require(raw['force_convention'] == 'negative derivative of the same native vacuum energy'
                and np.array_equal(item['MACE_forces_eV_A'], raw['forces_eV_A'])
                and np.allclose(composite, mace + solvent, rtol=1e-12, atol=1e-10)
                and np.allclose(-composite, item['composite_gradient']['kcal_mol_A'], rtol=1e-12, atol=1e-10),
                'Composite force sign/units or saved MACE vector differs')
        for name, force in zip(forces, (composite, mace, solvent)):
            forces[name][metal] = force
    for name, medium in (('GFN2_vacuum', 'vacuum'), ('GFN2_water', 'ALPB(Water)')):
        forces[name] = {}
        for metal in METALS:
            task = prepared['endpoints']['whole'][metal]
            gradient = np.asarray(raw_rows[(task['task_id'], medium)]['gradient']['hartree_per_bohr'])
            require(gradient.shape == (256, 3) and np.isfinite(gradient).all(), 'Invalid saved component gradient')
            forces[name][metal] = -gradient * gradient_units['hartree_to_kcal_mol'] / gradient_units['bohr_to_A']
    for metal in METALS:
        require(np.allclose(forces['GFN2_water'][metal] - forces['GFN2_vacuum'][metal], forces['solvent'][metal],
                            rtol=1e-12, atol=1e-10), 'Saved raw solvent gradient difference/units differs')
    for component in forces.values():
        component['common'] = (component['La'] + component['Dy']) / 2
        component['Dy_minus_La'] = component['Dy'] - component['La']
    loads = {c: {state: normalized.T @ f.reshape(-1) for state, f in states.items()} for c, states in forces.items()}
    rows = []
    for k, mode in enumerate(data['modes']):
        rows.append({'mode_index': k, 'mode_name': mode['name'], 'family': mode['family'], 'residue': mode.get('residue'),
            'coordinate_unit': mode['unit'], 'physical_column_norm': float(scales[k]),
            'attribution': 'backbone' if mode['family'] in ('phi', 'psi') else mode['family'],
            'normalized_force_loads_kcal_mol_A': {c: {state: float(v[k]) for state, v in states.items()} for c, states in loads.items()},
            'generalized_composite_force': {state: float(v[k] * scales[k]) for state, v in loads['composite'].items()},
            'generalized_force_units': 'kcal/mol/radian' if mode['unit'] == 'radian' else 'kcal/mol/angstrom'})
    rankings = {state: sorted(range(64), key=lambda k: (-abs(float(loads['composite'][state][k])), k))
                for state in ('La', 'Dy', 'common', 'Dy_minus_La')}
    families = {}
    for name, included in (('backbone', ('phi', 'psi')), ('sidechain', ('sidechain',)), ('metal_translation', ('metal_translation',))):
        indices = [k for k, m in enumerate(data['modes']) if m['family'] in included]
        families[name] = {'mode_indices': indices, 'summary_nonadditive': True,
            'composite_loads': {state: {'max_absolute_kcal_mol_A': float(np.max(np.abs(v[indices]))),
                                      'RMS_kcal_mol_A': float(np.sqrt(np.mean(v[indices] ** 2)))}
                                for state, v in loads['composite'].items()}}
    result = {'source_id': source, 'status': 'profiled', 'modes': mode_row['modes'],
        'source_preparation': prep_row['preparation'], 'physical_source_ids': ids,
        'force_units': 'kcal/mol/angstrom; positive is negative energy derivative', 'mode_rows': rows,
        'rankings_by_absolute_composite_load': rankings, 'family_attribution': families,
        'source_local_metric': metric(normalized, np),
        'normalization': 'equal-weight Euclidean norm over all256 physical atoms, including H; no cross-source Cartesian subtraction',
        'attribution_limit': 'overlapping coordinate loads and family RMS values are not additive energy or force partitions'}
    return dict(result=result, normalized=normalized, forces=forces, loads=loads, kin=kin,
                modes=modes, prepared=prepared, atoms=atoms, bonds=bonds, jac=jac, x=x)


def select_basis(by, np):
    modes = by[SOURCES[0]]['modes']['modes']
    require(by[SOURCES[0]]['modes']['mode_order'] == by[SOURCES[1]]['modes']['mode_order'], 'Source mode names/order differ')
    selected = [k for k, m in enumerate(modes) if m['family'] == 'metal_translation']
    require([modes[k]['name'] for k in selected] == RULE['mandatory'], 'Target translation catalogue differs')
    score = [max(abs(float(by[s]['loads']['composite'][state][k])) for s in SOURCES
                 for state in ('La', 'Dy', 'Dy_minus_La')) for k in range(64)]
    order = sorted((k for k, m in enumerate(modes) if m['unit'] == 'radian'), key=lambda k: (-score[k], k))
    ranking = []
    for k in order:
        residuals = {}
        for source in SOURCES:
            matrix = by[source]['normalized']; basis = matrix[:, selected]; vector = matrix[:, k]
            residuals[source] = float(np.linalg.norm(vector - basis @ np.linalg.lstsq(basis, vector, rcond=RCOND)[0]))
        accepted = len(selected) < 7 and all(r > RESIDUAL_MIN for r in residuals.values())
        row = {'mode_index': k, 'mode_name': modes[k]['name'], 'family': modes[k]['family'],
               'score_kcal_mol_A': score[k], 'residual_norm_by_source': residuals, 'selected': accepted,
               'reason': 'accepted' if accepted else ('four_torsion_limit' if len(selected) >= 7 else 'dependent_in_at_least_one_source')}
        ranking.append(row)
        if accepted:
            selected.append(k)
    for source in SOURCES:
        by[source]['result']['active_metric'] = metric(by[source]['normalized'][:, selected], np)
        by[source]['result']['active_own_and_differential_loads_kcal_mol_A'] = {
            state: values[selected].tolist() for state, values in by[source]['loads']['composite'].items()}
    return {'selection_rule': RULE, 'mode_indices_selection_order': selected,
            'mode_indices_kinematic_order': sorted(selected), 'mode_names_selection_order': [modes[k]['name'] for k in selected],
            'complete_torsion_ranking': ranking, 'common_names_different_source_local_Jacobians': True,
            'individual_metal_loads_retained': True, 'optimization_performed': False}


def finite_geometries(source, out, np, Geometry, geometry_settings):
    modes = source['modes']['modes']; jac = source['jac']; x = source['x']; kin = source['kin']
    differential = source['loads']['composite']['Dy_minus_La']
    k = min((k for k, m in enumerate(modes) if m['unit'] == 'radian'), key=lambda k: (-abs(float(differential[k])), k))
    mode = modes[k]; rmax = float(np.linalg.norm(jac[k], axis=1).max())
    require(rmax > max(STEPS_A) / 2, 'Degenerate finite-check torsion')
    atoms = source['atoms']; bonds = [b['native_particle_indices'] for b in source['bonds']]
    adapted = [dict(source_id=a['id'], residue_id=a['id'].rsplit('/', 1)[0], name=a['id'].rsplit('/', 1)[1],
                    element=a['element'], resname=a['resname']) for a in atoms[:255]]
    moving = mode['moving_indices']; require(255 not in moving, 'Derivative nominee moves target')
    guard = Geometry(adapted, bonds, x[:255], x[:255], moving, [{'element': 'La', 'xyz_A': x[255].tolist()}])
    origin_guard = guard.check(x[:255]); require(origin_guard['pass'], 'Actual source fails existing covalent/contact guards')
    qcheck = np.zeros(64); qcheck[k] = 2 * np.arcsin(min(STEPS_A) / (2 * rmax))
    derivative_check = kin.check(qcheck)
    require(derivative_check['pass'], 'Existing kinematic Jacobian/bond/fixed check failed')
    graph = [set() for _ in atoms]
    for a, b in bonds:
        graph[a].add(b); graph[b].add(a)
    angles = [(a, b, c) for b, ns in enumerate(graph) for a in sorted(ns) for c in sorted(ns) if a < c]
    def cosines(y):
        a, b, c = np.asarray(angles).T; u = y[a] - y[b]; v = y[c] - y[b]
        return np.einsum('ij,ij->i', u, v) / np.linalg.norm(u, axis=1) / np.linalg.norm(v, axis=1)
    origin_cosines = cosines(x)
    protected = [g['indices'] for g in source['modes']['protected_groups']]
    tangent = jac[k] / rmax
    result = {'source_id': '7CCO', 'mode_index': k, 'mode': mode, 'selection_rule': RULE['FD_nominee'],
        'normalized_differential_force_load_kcal_mol_A': float(differential[k]),
        'maximum_radius_A': rmax, 'step_parameter': 'signed maximum physical chord displacement h in angstrom',
        'q_formula': 'theta(h)=2*asin(h/(2*rmax)); all other coordinates zero',
        'tangent_dpositions_dh': tangent.tolist(),
        'central_difference_formula': '(E(+h)-E(-h))/(2*h); compare to dE/dh at h=0',
        'analytic_energy_derivatives_kcal_mol_A': {c: {state: float(-np.sum(f * tangent)) for state, f in states.items()}
                                                for c, states in source['forces'].items()},
        'numerical_admission_from_plan': {'absolute_kcal_mol_A': 0.2, 'relative': 0.05, 'steps_A': list(STEPS_A)},
        'independent_energy_derivative_check_performed': False, 'kinematics_check_at_positive_half_step': derivative_check,
        'source_guard': origin_guard, 'geometry_guard_settings': geometry_settings, 'geometries': [], 'tasks': []}
    for h in STEPS_A:
        for sign in (1, -1):
            name = ('h001' if h == 0.001 else 'h0005') + ('_plus' if sign > 0 else '_minus')
            directory = out / name; directory.mkdir(parents=True)
            q = np.zeros(64); q[k] = 2 * np.arcsin(sign * h / (2 * rmax)); y = kin.evaluate(q)[0]
            inherited = guard.check(y[:255]); displacement = float(np.linalg.norm(y - x, axis=1).max())
            angle_error = float(np.max(np.abs(cosines(y) - origin_cosines)))
            group_error = max(float(np.max(np.abs(np.linalg.norm(y[g][:, None] - y[g][None, :], axis=2) -
                                                  np.linalg.norm(x[g][:, None] - x[g][None, :], axis=2)))) for g in protected)
            fixed = sorted(set(range(256)) - set(moving))
            fixed_error = float(np.max(np.abs(y[fixed] - x[fixed])))
            passed = bool(inherited['pass'] and abs(displacement - h) <= 1e-12 and angle_error <= 1e-10
                          and group_error <= 1e-9 and fixed_error <= 1e-12
                          and inherited['source_bond_max_error_A'] <= 1e-9)
            record = {'geometry_id': name, 'status': 'prepared' if passed else 'unavailable', 'q': q.tolist(),
                      'mode_index': k, 'mode_name': mode['name'], 'theta_radian': float(q[k]),
                      'signed_step_A': sign * h, 'actual_maximum_physical_displacement_A': displacement,
                      'inherited_guard': inherited, 'max_bond_angle_cosine_change': angle_error,
                      'max_protected_group_distance_change_A': group_error, 'fixed_atom_max_change_A': fixed_error}
            physical = [dict(a, xyz_A=y[i].tolist()) for i, a in enumerate(atoms)]
            write(directory / 'PHYSICAL.json', physical); record['physical_source'] = pin(directory / 'PHYSICAL.json')
            write(directory / 'GEOMETRY.json', record)
            for metal in METALS:
                task = source['prepared']['endpoints']['whole'][metal]
                path = directory / (metal + '.xyz')
                with path.open('x') as handle:
                    handle.write('256\n7CCO whole source; paired exact torsion derivative check\n')
                    for i, atom in enumerate(atoms):
                        symbol = metal if atom['element'] == 'M' else atom['element']
                        handle.write(symbol + ' ' + ' '.join(format(v, '.17g') for v in y[i]) + '\n')
                result['tasks'].append(dict(task, task_id=task['task_id'] + '_FD_' + name, origin_task_id=task['task_id'],
                    origin_xyz=task['xyz'], origin_physical_source=task['physical_source'], xyz=pin(path),
                    physical_source=record['physical_source'], geometry=pin(directory / 'GEOMETRY.json'),
                    geometry_id=name, status=record['status'], signed_step_A=sign * h))
            result['geometries'].append(record)
    result['status'] = 'prepared' if all(g['status'] == 'prepared' for g in result['geometries']) else 'unavailable'
    result['no_geometry_or_energy_conditioned_retry'] = True
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('config', 'output', 'repo-root'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args(); require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'All numerical work requires Slurm')
    start = time.monotonic(); cfg = json.loads(args.config.read_text()); verify(cfg['plan'])
    if 'implementation' in cfg:
        require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Wrong adapter')
    modes = read(cfg['modes_result']); solvent = read(cfg['solvent_result'])
    require(modes['status'] == solvent['status'] == 'complete' and modes['preparation'] == solvent['preparation'],
            'Complete matched modes/solvent sources required')
    require(solvent['gradient_units']['eV_to_kcal_mol'] == EV_TO_KCAL and solvent['mace_forces'], 'Saved force units or MACE force receipt missing')
    prep = read(modes['preparation']); mace = read(solvent['mace_forces']); mcfg = read(mace['config'])
    require(prep['status'] == mace['status'] == 'complete' and mcfg['preparation'] == modes['preparation'], 'MACE force preparation differs')
    require(mace['versions'] == mcfg['expected_versions'] and mace['model'] == mcfg['model'], 'MACE model/runtime differs from admitted force call')
    corrections = unique(solvent['corrections'], lambda r: (r['source_id'], r['metal']))
    require(set(corrections) == {(s, m) for s in SOURCES for m in METALS}, 'Exact four composite force cells required')
    raw_rows = unique(solvent['rows'], lambda r: (r['task']['origin_task_id'], r['task']['medium']))
    require(len(raw_rows) == 8 and all(r['status'] == 'accepted' for r in raw_rows.values()), 'Eight accepted solvent derivative endpoints required')
    for item in corrections.values():
        for medium in ('vacuum', 'ALPB(Water)'):
            task = raw_rows[(item['origin_task_id'], medium)]['task']
            require(task['source_xyz'] == item['xyz'] and task['source_id'] == item['source_id'] and task['metal'] == item['metal']
                    and task['charge'] == -2 and task['multiplicity'] == 1 and task['uhf'] == 0
                    and task['physical_multiplicity'] == (1 if item['metal'] == 'La' else 6), 'GFN2 source/effective/physical state differs')
            verify(task['xyz'])
    mode_rows = unique(modes['rows'], lambda r: r['source_id']); prep_rows = unique(prep['rows'], lambda r: r['source_id'])
    force_by = unique(mace['rows'], lambda r: r['task_id'])
    require(set(mode_rows) == set(prep_rows) == set(SOURCES) and len(force_by) == 4, 'Source denominator differs')
    names = set(modes['helpers']) | {'scaffold_accommodation.py'}
    require(set(cfg['helpers']) == names, 'Exact kinematics and geometry helper pins required')
    for name in names:
        verify(cfg['helpers'][name])
        require(pin(args.repo_root / 'scripts' / name)['sha256'] == cfg['helpers'][name]['sha256'], 'Current helper differs: ' + name)
        if name in modes['helpers']:
            require(cfg['helpers'][name]['sha256'] == modes['helpers'][name]['sha256'], 'Helper differs from saved modes: ' + name)
    sys.path.insert(0, str(args.repo_root.resolve() / 'scripts'))
    import numpy as np
    import scipy
    from mace_site_kinematics import Kinematics
    from scaffold_accommodation import Geometry, SETTINGS
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    result = {'protocol_id': 'lbt3_saved_composite_force_profile_v1', 'status': 'incomplete',
        'config': pin(args.config), 'plan': cfg['plan'], 'implementation': pin(__file__), 'helpers': cfg['helpers'],
        'modes_result': cfg['modes_result'], 'solvent_result': cfg['solvent_result'], 'mace_forces': solvent['mace_forces'],
        'selection_rule': RULE, 'new_molecular_calls': 0, 'optimizer_steps': 0, 'energy_scoring_performed': False,
        'versions': {'numpy': np.__version__, 'scipy': scipy.__version__}, 'rows': [], 'slurm_job_id': os.environ['SLURM_JOB_ID']}
    try:
        by = {s: source_profile(s, mode_rows[s], prep_rows[s], corrections, force_by, raw_rows,
                                solvent['gradient_units'], np, Kinematics) for s in SOURCES}
        result['active_candidate_basis'] = select_basis(by, np)
        result['rows'] = [by[s]['result'] for s in SOURCES]
        write(out / 'PROFILE.json', dict(result, status='profiled'))
        result['profile'] = pin(out / 'PROFILE.json')
        finite = finite_geometries(by['7CCO'], out / 'finite_check', np, Geometry, SETTINGS)
        write(out / 'FINITE_CHECK.json', finite); result['finite_check'] = pin(out / 'FINITE_CHECK.json')
        result['finite_geometries_generated'] = len(finite['geometries'])
        result['paired_finite_tasks_prepared'] = sum(t['status'] == 'prepared' for t in finite['tasks'])
        result['status'] = 'complete' if finite['status'] == 'prepared' else 'geometry_unavailable'
    except Exception as exc:
        result.update(status='unavailable', failure=str(exc), exception=type(exc).__name__, traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic() - start
    result['allocation'] = {k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_CPU')}
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}), flush=True)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
