"""One own-gradient eight-coordinate proposal per Hans source/metal, on Slurm.

Reuse saved gradients and the complete existing classical ledger. No quantum
calls, line search, simulation, new parameterization, or differential optimizer.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import socket
import sys
import time
import traceback

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
SOURCES = ('Hans8DQ2', 'Hans8FNR')
METALS = ('La', 'Dy')
SETTINGS = dict(max_atom_displacement_A=.02, max_torsion_degrees=.5,
                gram_relative_cutoff=1e-10, max_bond_change_A=1e-10,
                max_angle_change_rad=1e-10, max_signed_volume_change_A3=1e-9)
# Ancestors precede their descendants; the derivatives of the composed tree
# rotations are the current-axis cross products on their complete distal sets.
AXES = (('N83_CA83', 'A/83//N', 'A/83//CA'),
        ('Asp85_CA_CB', 'A/85//CA', 'A/85//CB'),
        ('Asp85_CB_CG', 'A/85//CB', 'A/85//CG'),
        ('Asp87_CA_CB', 'A/87//CA', 'A/87//CB'),
        ('Asp87_CB_CG', 'A/87//CB', 'A/87//CG'))


def pin(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def write_new(path, value):
    text = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n'
    with Path(path).open('x') as handle:
        handle.write(text)


def prepare(config_path, output):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Geometry and classical preparation require Slurm')
    import numpy as np
    import openmm as mm
    from openmm import unit
    from scipy.spatial import cKDTree
    sys.path.insert(0, str(ROOT / 'scripts'))
    from affordable_common import read_json, verify, InvalidArtifact, xyz
    from affordable_response import cap_jacobians
    from metal_environment_force_assembly import prepare_mapping, target_source_id
    from metal_environment_component_checks import build_system, GROUPS
    from metal_environment_lanm_reference import check_config
    from scaffold_accommodation import Geometry, SETTINGS as GEOMETRY_SETTINGS

    dependencies = [pin(m.__file__) for m in list(sys.modules.values())
                    if getattr(m, '__file__', None) and Path(m.__file__).resolve().parent == ROOT / 'scripts']
    consumed = {}

    def require(condition, message):
        if not condition:
            raise InvalidArtifact(message)

    def read_pin(record):
        path = verify(record)
        consumed[str(path.resolve())] = pin(path)
        return read_json(path)

    def unique(atoms):
        result = {a['id']: a for a in atoms}
        require(len(result) == len(atoms), 'Duplicate physical IDs')
        return result

    cfg = read_json(config_path)
    read_pin(cfg['profile'])
    verify(cfg['plan'])
    consumed[cfg['plan']['path']] = cfg['plan']
    profile = read_pin(cfg['profile'])
    final = read_pin(profile['final_comparison'])
    final_rows = {(r['source'], r['metal']): r for r in final['rows']}
    result = dict(protocol='lanm_coupled_eight_coordinate_preparation_v1', config=pin(config_path),
        settings=SETTINGS, mode_order=[r[0] for r in AXES] + ['target_x', 'target_y', 'target_z'],
        mode_units=['radian'] * 5 + ['angstrom'] * 3, sources={},
        new_quantum_calls=0, integration_steps=0, actual_component_queries=0,
        classifier=None, affinity=None, full_hybrid_qualified=False,
        common_pool_policy='origin plus both own-metal proposals; paired evaluation, no differential optimization',
        limitations=['Fixed monomer source states and unequal crystal water inventories remain',
                     'Whole native finite classical ledger retained; bulk solvent and dimer absent',
                     'Linear electronic work is a diagnostic, never a score or finite quantum energy',
                     'Full hybrid forces and cross-LJ transfer remain unqualified'])

    for source in SOURCES:
        source_out = output / source
        source_out.mkdir()
        assembled = {m: read_pin(cfg['assemblies'][source][m]) for m in METALS}
        base = assembled['La']
        inputs_pin = base['inputs']
        inputs = read_pin(inputs_pin)
        export = read_pin(inputs['export'])
        native = read_pin(export['atoms'])
        original_physical = read_pin(inputs['physical_sources']['A'])
        bonds = read_pin(export['bonds'])
        core = read_pin(inputs['configurations']['A']['core_mapping'])
        field = read_pin(inputs['environment_atoms'])
        boundary = read_pin(inputs['boundary_mapping'])
        target_ids = [target_source_id(a) for a in core if a['kind'] == 'metal']
        require(len(target_ids) == 1, 'Expected exactly one EF3 target')
        target = target_ids[0]
        atoms = base['physical_atoms']
        ids = [a['id'] for a in atoms]
        ix = {sid: i for i, sid in enumerate(ids)}
        require(len(ix) == len(ids), 'Duplicate assembled IDs')
        x = np.asarray([a['xyz_A'] for a in atoms])
        metadata = unique(native)
        protein_ids = [a['id'] for a in native if a['kind'] == 'protein_source']
        protein_set = set(protein_ids)
        require(len(protein_ids) == 1665, 'Hans protein inventory changed')
        gradient = {}
        old_classical = {}
        for metal in METALS:
            a = assembled[metal]
            f = final_rows[(source, metal)]
            require(a['inputs'] == inputs_pin == f['inputs'] and a['collection'] == f['collection']
                    and a['classical_record'] == f['classical'] and a['metal'] == metal,
                    'Assembled origin differs from accepted six-cell provenance')
            require([v['id'] for v in a['physical_atoms']] == ids and
                    np.array_equal([v['xyz_A'] for v in a['physical_atoms']], x),
                    'Paired origin physical ordering/coordinates differ')
            gradient[metal] = np.asarray(a['total_gradient_kcal_mol_A'])
            require(gradient[metal].shape == x.shape and np.isfinite(gradient[metal]).all(), 'Invalid total gradient')
            require(np.array_equal(-gradient[metal], a['total_force_kcal_mol_A']), 'Gradient/force sign mismatch')
            old_classical[metal] = read_pin(a['classical_record'])
        gradients = dict(gradient, common=(gradient['La'] + gradient['Dy']) / 2,
                         Dy_minus_La=gradient['Dy'] - gradient['La'])
        ledger_pin = old_classical['La']['ledger']
        require(ledger_pin == old_classical['Dy']['ledger'], 'Paired origins use different classical ledgers')
        ledger = read_pin(ledger_pin)
        require(ledger['inputs'] == inputs_pin and ledger['native_system'] == export['system'],
                'Classical ledger does not bind this native source')
        require(ledger['boundary_policy'] == 'fixed_cap_radial_bonds_retained_substitutable_angular_terms_omitted_v1',
                'Boundary term policy changed')
        artifacts = {name: read_pin(rec) for name, rec in ledger['artifacts'].items()}
        ledger_ids = [target if a['id'] == 'target_EF3' else a['id'] for a in artifacts['particles.json']]
        require(len(ledger_ids) == len(set(ledger_ids)) and set(ledger_ids) == set(ids), 'Ledger inventory differs')
        ledger_order = [ix[sid] for sid in ledger_ids]
        require(np.allclose([a['xyz_A'] for a in artifacts['particles.json']], x[ledger_order], atol=1e-8, rtol=0),
                'Ledger origin coordinates differ from saved full gradients')
        require(pin(ROOT / 'scripts/metal_environment_component_checks.py') == ledger['reused_builder'],
                'Native classical builder changed from accepted ledger')
        graph = {sid: set() for sid in ids}
        for a, b in bonds:
            require(a in ix and b in ix, 'Native bond lacks physical mapping')
            graph[a].add(b)
            graph[b].add(a)
        modes = []
        for name, a, b in AXES:
            require(a in graph and b in graph[a], 'Mode axis is not a native source bond')
            seen, pending = {b}, [b]
            while pending:
                current = pending.pop()
                for neighbor in graph[current] - seen:
                    if {current, neighbor} == {a, b}:
                        continue
                    seen.add(neighbor)
                    pending.append(neighbor)
            require(a not in seen and seen <= protein_set, 'Torsion does not cut a protein graph bridge')
            modes.append(dict(name=name, axis=[a, b], moving_source_ids=sorted(seen),
                              moving_indices=[ix[sid] for sid in sorted(seen)]))
        for i, parent in enumerate(modes):
            pset = set(parent['moving_source_ids'])
            for child in modes[i+1:]:
                cset = set(child['moving_source_ids'])
                require(not pset & cset or cset <= pset and set(child['axis']) <= pset,
                        'Mode order is not ancestor-before-descendant/disjoint')

        def jacobian(y):
            j = np.zeros((len(ids), 3, 8))
            for k, mode in enumerate(modes):
                a, b = [ix[sid] for sid in mode['axis']]
                axis = y[b] - y[a]
                axis /= np.linalg.norm(axis)
                moving = mode['moving_indices']
                j[moving, :, k] = np.cross(axis, y[moving] - y[b])
            j[ix[target], :, 5:] = np.eye(3)
            return j

        def coordinates(q):
            y = x.copy()
            for k, mode in enumerate(modes):
                a, b = [ix[sid] for sid in mode['axis']]
                axis = y[b] - y[a]
                axis /= np.linalg.norm(axis)
                moving = mode['moving_indices']
                pivot = y[b].copy()
                v = y[moving] - pivot
                y[moving] = pivot + v*np.cos(q[k]) + np.cross(axis, v)*np.sin(q[k]) + np.outer(v@axis, axis)*(1-np.cos(q[k]))
            y[ix[target]] += q[5:]
            return y

        j0 = jacobian(x)
        matrix = j0.reshape(-1, 8)
        scales = np.linalg.norm(matrix, axis=0)
        require(np.all(scales > 0), 'Zero coordinate Jacobian column')
        normalized = matrix / scales
        gram = normalized.T @ normalized
        eigenvalues, eigenvectors = np.linalg.eigh(gram)
        retained = eigenvalues > SETTINGS['gram_relative_cutoff'] * eigenvalues[-1]
        require(np.any(retained), 'No retained physical coordinate span')
        inverse = (eigenvectors[:, retained] / eigenvalues[retained]) @ eigenvectors[:, retained].T
        origin_loads = {name: (normalized.T @ g.reshape(-1)).tolist() for name, g in gradients.items()}
        mode_record = dict(physical_source_ids=ids, axes=modes, target_source_id=target,
            mode_order=result['mode_order'], mode_units=result['mode_units'],
            origin_jacobian=j0.tolist(), column_norms=scales.tolist(), normalized_Gram=gram.tolist(),
            Gram_eigenvalues=eigenvalues.tolist(), retained_span_rank=int(retained.sum()),
            retained_Gram_condition=float(eigenvalues[retained].max()/eigenvalues[retained].min()),
            relative_cutoff=SETTINGS['gram_relative_cutoff'], normalized_generalized_gradient_loads=origin_loads,
            derivative_formula='Each graph torsion: current unit axis cross (atom-current pivot); target translations: identity',
            whole_classical_ledger=ledger_pin)
        write_new(source_out / 'MODES.json', mode_record)

        # Reuse source chirality/peptide and relative severe-contact gates. No
        # parent energy model is constructed by this geometry-only helper.
        protein_index = {sid: i for i, sid in enumerate(protein_ids)}
        protein_atoms = [dict(source_id=sid, residue_id=sid.rsplit('/', 1)[0],
            name=sid.rsplit('/', 1)[1], element=metadata[sid]['element'], resname=metadata[sid]['resname'])
            for sid in protein_ids]
        protein_bonds = [(protein_index[a], protein_index[b]) for a, b in bonds if a in protein_set and b in protein_set]
        protein_order = [ix[sid] for sid in protein_ids]
        moving_ids = set().union(*(set(m['moving_source_ids']) for m in modes)) | {target}
        moving_indices = {ix[sid] for sid in moving_ids}
        free_protein = [protein_index[sid] for sid in protein_ids if sid in moving_ids]
        other = [dict(a, xyz_A=x[ix[a['id']]].tolist()) for a in atoms if a['id'] not in protein_set]
        guard = Geometry(protein_atoms, protein_bonds, x[protein_order], x[protein_order], free_protein, other)
        all_bonds = [(ix[a], ix[b]) for a, b in bonds]
        all_angles = [(ix[a], ix[center], ix[b]) for center, neighbors in graph.items()
                      for a in sorted(neighbors) for b in sorted(neighbors) if a < b]
        excluded = set()
        for sid, neighbors in graph.items():
            near = set(neighbors)
            for neighbor in neighbors:
                near.update(graph[neighbor])
            excluded.update(tuple(sorted((ix[sid], ix[n]))) for n in near if n != sid)

        def angle(y, a, b, c):
            u, v = y[a] - y[b], y[c] - y[b]
            return np.arccos(np.clip(u@v/np.linalg.norm(u)/np.linalg.norm(v), -1, 1))

        def check_geometry(y, q):
            guard.other = np.asarray([y[ix[a['id']]] for a in other])
            inherited = guard.check(y[protein_order])
            bond_error = max(abs(np.linalg.norm(y[a]-y[b])-np.linalg.norm(x[a]-x[b])) for a, b in all_bonds)
            angle_error = max(abs(angle(y, a, b, c)-angle(x, a, b, c)) for a, b, c in all_angles)
            volume_error, inversions = 0., []
            for sid, neighbors in graph.items():
                if len(neighbors) != 4:
                    continue
                order = [ix[n] for n in sorted(neighbors)]
                p0, p1 = x[order], y[order]
                v0, v1 = np.linalg.det(p0[:3]-p0[3]), np.linalg.det(p1[:3]-p1[3])
                volume_error = max(volume_error, abs(float(v1-v0)))
                if v0*v1 < 0:
                    inversions.append(sid)
            contacts = []
            for a, b in sorted(cKDTree(y).query_pairs(GEOMETRY_SETTINGS['new_overlap_heavy_A'])):
                if (a not in moving_indices and b not in moving_indices) or (a, b) in excluded:
                    continue
                cutoff = GEOMETRY_SETTINGS['new_overlap_H_A'] if 'H' in (atoms[a]['element'], atoms[b]['element']) else GEOMETRY_SETTINGS['new_overlap_heavy_A']
                d, d0 = float(np.linalg.norm(y[a]-y[b])), float(np.linalg.norm(x[a]-x[b]))
                if d < min(cutoff, d0)-GEOMETRY_SETTINGS['overlap_roundoff_A']:
                    contacts.append(dict(source_ids=[ids[a], ids[b]], source_A=d0, candidate_A=d, cutoff_A=cutoff))
            fixed = [i for i in range(len(ids)) if i not in moving_indices]
            fixed_error = float(np.max(np.abs(y[fixed]-x[fixed])))
            max_move = float(np.linalg.norm(y-x, axis=1).max())
            max_angle = float(np.rad2deg(np.abs(q[:5])).max())
            admitted = bool(inherited['pass'] and bond_error <= SETTINGS['max_bond_change_A']
                and angle_error <= SETTINGS['max_angle_change_rad'] and volume_error <= SETTINGS['max_signed_volume_change_A3']
                and not inversions and not contacts and fixed_error <= 1e-12
                and max_move <= SETTINGS['max_atom_displacement_A']+1e-12
                and max_angle <= SETTINGS['max_torsion_degrees']+1e-12)
            return dict(admitted=admitted, inherited_geometry=inherited, max_bond_change_A=bond_error,
                max_angle_change_rad=angle_error, max_signed_volume_change_A3=volume_error,
                signed_volume_inversions=inversions, worsened_severe_contacts=contacts,
                fixed_atom_max_error_A=fixed_error, max_atom_displacement_A=max_move, max_torsion_degrees=max_angle,
                geometry_rule_source=pin(ROOT/'scripts/scaffold_accommodation.py'),
                contact_rule='reject distance < min(existing cutoff, same-pair origin distance) - 1e-7 A; exclude 1-2/1-3')

        def write_candidate(directory, y, q, proposal_metal):
            current_physical = copy.deepcopy(original_physical)
            for a in current_physical:
                a['xyz_A'] = y[ix[a['id']]].tolist()
            write_new(directory/'physical_source.json', current_physical)
            current_native = copy.deepcopy(native)
            for a in current_native:
                a['xyz_A'] = y[ix[a['id']]].tolist()
            write_new(directory/'atoms_charges.json', current_native)
            current_export = copy.deepcopy(export)
            current_export.update(atoms=pin(directory/'atoms_charges.json'), origin_export=inputs['export'],
                coordinate_changes='declared eight-coordinate candidate; native topology/charges unchanged')
            write_new(directory/'EXPORT.json', current_export)
            env = copy.deepcopy(field)
            for a in env:
                a['xyz_A'] = y[ix[a['id']]].tolist()
            write_new(directory/'environment_atoms.json', env)
            with (directory/'environment.pc').open('x') as handle:
                handle.write(str(len(env))+'\n')
                handle.write(''.join(' '.join(format(z, '.17g') for z in (a['charge_e'], *a['xyz_A']))+'\n' for a in env))
            current_boundary = copy.deepcopy(boundary)
            current_boundary['original_boundary_mapping'] = inputs['boundary_mapping']
            current_boundary['coordinate_preparation'] = 'own_total_gradient_coupled_candidate_v1'
            for sid in current_boundary['retained_physical_CA_coordinates']:
                current_boundary['retained_physical_CA_coordinates'][sid] = y[ix[sid]].tolist()
            for row in current_boundary['ledger']:
                row['original_dipole_diagnostic'] = {k: row[k] for k in ('dipole_before_eA', 'retained_MM_dipole_after_eA')}
                prefix = row['partial_residue']
                row['dipole_before_eA'] = sum((a['charge_e']*np.asarray(a['xyz_A']) for a in current_native if a['id'].startswith(prefix)), np.zeros(3)).tolist()
                row['retained_MM_dipole_after_eA'] = sum((a['charge_e']*np.asarray(a['xyz_A']) for a in env if a['id'].startswith(prefix)), np.zeros(3)).tolist()
            write_new(directory/'boundary_mapping.json', current_boundary)
            mapping = copy.deepcopy(core)
            j = jacobian(y)
            core_j = np.zeros((len(core), 3, 8))
            for a in mapping:
                qi = a['qm_index']
                if a['kind'] == 'metal':
                    a['xyz_A'] = y[ix[target]].tolist()
                    a['source_site']['xyz_A'] = a['xyz_A']
                    core_j[qi] = j[ix[target]]
                elif a['kind'] == 'source':
                    a['xyz_A'] = y[ix[a['source_id']]].tolist()
                    core_j[qi] = j[ix[a['source_id']]]
                else:
                    ri, oi = ix[a['retained_source_id']], ix[a['omitted_source_id']]
                    r, o = y[ri], y[oi]
                    jr, jo = cap_jacobians(r, o, a['length_A'])
                    a.update(xyz_A=(r+a['length_A']*(o-r)/np.linalg.norm(o-r)).tolist(),
                        retained_xyz_A=r.tolist(), omitted_xyz_A=o.tolist(),
                        jacobian_retained=jr.tolist(), jacobian_omitted=jo.tolist())
                    core_j[qi] = jr@j[ri] + jo@j[oi]
            write_new(directory/'core_mapping.json', mapping)
            endpoints = {}
            for metal in METALS:
                endpoint = copy.deepcopy(inputs['configurations']['A']['endpoints'][metal])
                old_xyz = xyz(verify(endpoint['xyz']))
                path = directory/(metal+'.xyz')
                with path.open('x') as handle:
                    handle.write(str(len(mapping))+'\n'+source+' '+proposal_metal+' own-gradient candidate\n')
                    handle.write(''.join(a[0]+' '+' '.join(format(z, '.17g') for z in m['xyz_A'])+'\n' for a, m in zip(old_xyz, mapping)))
                endpoint['xyz'] = pin(path)
                endpoints[metal] = endpoint
            configuration = dict(endpoints=endpoints, core_mapping=pin(directory/'core_mapping.json'),
                                 pointcharges=pin(directory/'environment.pc'), geometry_id=source+'_'+proposal_metal+'_step')
            current_inputs = copy.deepcopy(inputs)
            current_inputs.update(protocol_id=inputs['protocol_id']+'_coupled_singleton_v1',
                origin_inputs=inputs_pin, origin_configuration='A', export=pin(directory/'EXPORT.json'),
                physical_source=pin(directory/'physical_source.json'),
                physical_sources={k: pin(directory/'physical_source.json') for k in ('A','B')},
                environment_atoms=pin(directory/'environment_atoms.json'), boundary_mapping=pin(directory/'boundary_mapping.json'),
                core_mapping={k: pin(directory/'core_mapping.json') for k in ('A','B')},
                configurations={k: copy.deepcopy(configuration) for k in ('A','B')},
                coordinate_preparation='own_total_gradient_coupled_candidate_v1',
                perturbation=dict(type='eight-coordinate coupled singleton', mode_order=result['mode_order'],
                    coordinates=q.tolist(), coordinate_units=result['mode_units'], selector_metal=proposal_metal,
                    selector='own total physical gradient only; no energy or label selection'),
                configuration_aliases={'B':'A'},
                alias_reason='Both labels denote one candidate for existing singleton A staging; do not calculate duplicate B',
                implementation=pin(__file__), preparation_plan=cfg['plan'], energy_force_calls=0,
                full_hybrid_status='matching classical candidate available separately; new electronic energy unavailable')
            write_new(directory/'INPUTS.json', current_inputs)
            check_config(current_inputs)
            for metal in METALS:
                mapped = prepare_mapping(directory/'INPUTS.json', 'A', metal)
                actual = unique(mapped['physical_atoms'])
                require(set(actual) == set(ids), 'Candidate physical mapping lost atoms')
                require(all(np.allclose(actual[sid]['xyz_A'], y[ix[sid]], atol=1e-8, rtol=0) for sid in ids),
                        'Candidate physical mapping differs, including zero-charge boundary anchors')
            write_new(directory/'KINEMATICS.json', dict(physical_source_ids=ids, coordinates=q.tolist(),
                mode_order=result['mode_order'], mode_units=result['mode_units'], physical_jacobian=j.tolist(),
                core_jacobian=core_j.tolist(), field_jacobian=j[[ix[a['id']] for a in env]].tolist(),
                cap_chain_rule='J_cap = J_retained*dR_retained/dq + J_omitted*dR_omitted/dq',
                no_synthetic_cap_degree_of_freedom=True))
            return pin(directory/'INPUTS.json')

        proposals = {}
        for metal in METALS:
            directory = source_out/(metal+'_step')
            directory.mkdir()
            coefficients = -inverse @ (normalized.T @ gradient[metal].reshape(-1))
            qdirection = coefficients/scales
            # Ancestor rotations preserve child-axis radii. Sum of arc-length
            # bounds is conservative for the composed finite displacement and
            # gives one scale without backtracking or post-outcome tuning.
            atom_upper = np.sum(np.linalg.norm(j0, axis=1)*np.abs(qdirection), axis=1)
            require(float(atom_upper.max()) > 0, 'Own gradient has zero retained-span displacement')
            bound_scale = SETTINGS['max_atom_displacement_A']/float(atom_upper.max())
            torsion_max = float(np.abs(qdirection[:5]).max())
            if torsion_max > 0:
                bound_scale = min(bound_scale, float(np.deg2rad(SETTINGS['max_torsion_degrees']))/torsion_max)
            q = bound_scale*qdirection
            tangent_step = (matrix@q).reshape(x.shape)
            require(float(np.sum(gradient[metal]*tangent_step)) < 0, 'Own total projected step is not descent')
            y = coordinates(q)
            geometry = check_geometry(y, q)
            write_new(directory/'PHYSICAL.json', dict(source_ids=ids, xyz_A=y.tolist(), origin=inputs_pin))
            write_new(directory/'GEOMETRY.json', geometry)
            proposal = dict(status='prepared' if geometry['admitted'] else 'unavailable_geometry',
                selector_metal=metal, q=q.tolist(), scale=bound_scale,
                single_scale_rule='minimum of Cartesian sum-of-arc upper-bound limit and per-torsion angle limit',
                geometry=pin(directory/'GEOMETRY.json'), physical=pin(directory/'PHYSICAL.json'),
                gradient_projection_work_diagnostic_kcal_mol={name: float(np.sum(g*tangent_step)) for name,g in gradients.items()},
                own_electronic_tangent_work_diagnostic_kcal_mol=float(np.sum(np.asarray(assembled[metal]['electronic']['gradient_kcal_mol_A'])*tangent_step)),
                own_classical_tangent_work_diagnostic_kcal_mol=float(np.sum(np.asarray(assembled[metal]['classical_gradient_kcal_mol_A'])*tangent_step)),
                finite_quantum_energy=None, finite_total_energy=None, quantum_ready=False,
                no_retry=True, classical={})
            if geometry['admitted']:
                proposal['inputs'] = write_candidate(directory, y, q, metal)
                proposal['kinematics'] = pin(directory/'KINEMATICS.json')
            proposals[metal] = proposal

        # Full parent-ledger mechanics for every candidate/target pair. Contexts
        # never integrate, and the original four classical origins are reused.
        for evaluated_metal in METALS:
            system, description = build_system(ledger, artifacts, evaluated_metal)
            integrator = mm.VerletIntegrator(.001)
            context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
            for proposer, proposal in proposals.items():
                if proposal['status'] != 'prepared':
                    continue
                directory = source_out/(proposer+'_step')
                physical = read_pin(proposal['physical'])
                y = np.asarray(physical['xyz_A'])
                context.setPositions(y[ledger_order]*.1*unit.nanometer)
                energies, component_gradients = {}, {}
                for group, name in enumerate(GROUPS):
                    state = context.getState(getEnergy=True, getForces=True, groups={group})
                    result['actual_component_queries'] += 1
                    energies[name] = float(state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
                    component_gradients[name] = -np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom))
                total = sum(component_gradients.values())
                require(np.isfinite(total).all() and np.isfinite(list(energies.values())).all(), 'Nonfinite classical candidate')
                path = directory/(evaluated_metal+'_CLASSICAL.json')
                record = dict(source_ids=ledger_ids, coordinates_A=y[ledger_order].tolist(), quantity='gradient',
                    units='kcal/mol/angstrom', metal=evaluated_metal, inputs=proposal['inputs'], origin_inputs=inputs_pin,
                    ledger=ledger_pin, target_source_id=target, QM_MM_Coulomb_included=False, C4_induction_included=False,
                    components_kcal_mol=energies, component_gradients_kcal_mol_A={k:v.tolist() for k,v in component_gradients.items()},
                    gradient_kcal_mol_A=total.tolist(), total_kcal_mol=float(sum(energies.values())),
                    builder_description=description, origin_classical=assembled[evaluated_metal]['classical_record'],
                    component_changes_kcal_mol={k:energies[k]-old_classical[evaluated_metal]['components_kcal_mol'][k] for k in GROUPS},
                    actual_classical_change_kcal_mol=float(sum(energies.values())-old_classical[evaluated_metal]['total_kcal_mol']),
                    finite_quantum_energy=None, finite_total_energy=None)
                write_new(path, record)
                proposal['classical'][evaluated_metal] = pin(path)
            del context, integrator
        write_new(source_out/'PROPOSALS.json', proposals)
        result['sources'][source] = dict(modes=pin(source_out/'MODES.json'), proposals=pin(source_out/'PROPOSALS.json'),
            accepted_geometry_count=sum(p['status']=='prepared' for p in proposals.values()),
            new_classical_cells=sum(len(p['classical']) for p in proposals.values()), reused_origin_classical_cells=2,
            harness_limitations=['Use explicit singleton A; B is an alias, not a new source or molecular cell',
                'Existing single-torsion collector must be replaced or explicitly adapted by root',
                'Moved retained_physical_CA_coordinates require explicit geometric seed handling'])
    for record in [*consumed.values(), *dependencies]:
        verify(record)
    result['consumed_pins'] = list(consumed.values())
    result['runtime_dependencies'] = dependencies
    result['runtime_versions'] = dict(python=sys.version, numpy=np.__version__, openmm=mm.__version__)
    write_new(output/'RESULT.json', result)
    return dict(result=pin(output/'RESULT.json'), actual_component_queries=result['actual_component_queries'],
                new_quantum_calls=0, integration_steps=0)


def main(args):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Numerical preparation requires Slurm')
    if int(os.environ.get('SLURM_NTASKS','1')) != 1 or int(os.environ.get('SLURM_CPUS_PER_TASK','1')) != 1:
        raise RuntimeError('Expected one task and one CPU')
    args.output.mkdir(parents=True, exist_ok=False)
    receipt = dict(status='running', config=pin(args.config), implementation=pin(__file__),
        slurm_job_id=os.environ['SLURM_JOB_ID'], host=socket.gethostname(), started_epoch=time.time(),
        new_quantum_calls=0, integration_steps=0, allocated_CPUs=1)
    write_new(args.output/'START.json', receipt)
    code = 1
    try:
        receipt.update(prepare(args.config, args.output))
        for key in ('config','implementation'):
            if pin(receipt[key]['path']) != receipt[key]:
                raise RuntimeError('Changed '+key+' during preparation')
        receipt['status'] = 'complete'
        code = 0
    except Exception as error:
        receipt.update(status='failed', error=str(error), traceback=traceback.format_exc(), result_accepted=False)
    finally:
        receipt['finished_epoch'] = time.time()
        receipt['elapsed_seconds'] = receipt['finished_epoch']-receipt['started_epoch']
        receipt['scheduler_accounting'] = 'collect separately after termination'
        write_new(args.output/'RECEIPT.json', receipt)
    return code


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.config, args.output = args.config.resolve(), args.output.resolve()
    raise SystemExit(main(args))
