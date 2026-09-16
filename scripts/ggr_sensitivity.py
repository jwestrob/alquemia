"""Approved GGR Stage C: frozen physical directions and analytic-gradient checks.

Prepares only the named center/directional tasks; never executes ORCA. Numerical
relaxation and entropy corrections remain unavailable. All numerical settings
and conditional half-step selection come from the pinned approved plan.
"""
from __future__ import annotations

import argparse
import copy
import math
from pathlib import Path
import re
import shutil

import gemmi
import numpy as np

import carve_generic as generic
from affordable_common import (BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, cache_key,
                               digest, energy, paired, read_json, record, verify,
                               write_new, xyz)
from affordable_peptide import SourceGraph
from affordable_response import cap_jacobians, extract, source_key

PROTOCOL = 'ggr_local_sensitivity_tightscf_dev_v1'
ROOT = Path(__file__).resolve().parents[1]
ORCA = Path('/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca')
REPRESENTATIONS = {
    'extended': 'generic_peptide_alpha_caps_native_r2scan3c_dev_v1',
    'connected': 'ggr_connected_segment_native_r2scan3c_dev_v1',
}
AMPLITUDES = {'metal': .02, 'peptide': math.pi / 180}
METHOD = 'r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF'
UNAVAILABLE = {'response_status': 'response_model_not_validated',
               'relaxation_correction_kcal_mol': None, 'entropy_correction_kcal_mol': None}


def input_state(path, *, gradient=None, tight=None):
    """Accept the exact native two-line scientific recipe, with no hidden blocks."""
    lines = [line.strip() for line in Path(path).read_text().splitlines()
             if line.strip() and not line.lstrip().startswith('#')]
    if len(lines) != 2 or not lines[0].startswith('!'):
        raise InvalidArtifact(f'unsupported scientific input recipe: {path}')
    tokens = lines[0][1:].lower().split()
    expected = {'r2scan-3c', 'noautostart', 'cpcm(water)', 'defgrid3'}
    if not expected.issubset(tokens) or set(tokens) - expected - {'tightscf', 'engrad'} or len(tokens) != len(set(tokens)):
        raise InvalidArtifact('native r2SCAN-3c/CPCM/DefGrid3 recipe required')
    for keyword, setting in [('engrad', gradient), ('tightscf', tight)]:
        if setting is not None and (keyword in tokens) != setting:
            raise InvalidArtifact(f'wrong {keyword} task type')
    declaration = re.fullmatch(r'\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+(\S+)', lines[1], re.I)
    if not declaration or int(declaration[2]) != 1:
        raise InvalidArtifact('singlet XYZ declaration required')
    return {'charge': int(declaration[1]), 'multiplicity': int(declaration[2]),
            'xyz_name': declaration[3], 'gradient': 'engrad' in tokens,
            'tightscf': 'tightscf' in tokens}


def atom_key(meta):
    return tuple(meta[k] for k in ('chain_index', 'residue_index', 'atom'))


def rotation(vector, axis, angle):
    v, a = np.asarray(vector), np.asarray(axis)
    return v * math.cos(angle) + np.cross(a, v) * math.sin(angle) + a * np.dot(a, v) * (1 - math.cos(angle))


def physical_model(repair):
    """Recover the actual bonded Gln140--Ile141 unit and source atom measure."""
    graph = SourceGraph(verify(repair['source_structure']), verify(repair['topology_definition']))
    carbonyl = graph.locate({'chain': 'A', 'resnum': 140, 'resname': 'GLN'})
    c, o, ca = [graph.key(carbonyl.key, name) for name in ('C', 'O', 'CA')]
    n = graph.amide_next.get(carbonyl.key)
    if n is None or graph.meta[n]['resname'] != 'ILE' or graph.meta[n]['resnum'] != 141:
        raise InvalidArtifact('approved Gln140--Ile141 peptide is absent')
    nb = graph.key(n[:2], 'CA')
    hydrogens = [h for h, parent in graph.hparents.items() if parent == n]
    if len(hydrogens) != 1:
        raise InvalidArtifact('approved amide-N hydrogen missing/ambiguous')
    moving = [c, o, n, hydrogens[0]]
    mapping = repair['atom_graph']['source_to_qm']
    selected = {atom_key(row['source']) for row in mapping if row['kind'] == 'source'}
    if not set(moving + [ca, nb]).issubset(selected):
        raise InvalidArtifact('physical crankshaft requires both retained alpha anchors and complete peptide')
    if len(selected) != sum(row['kind'] == 'source' for row in mapping):
        raise InvalidArtifact('duplicate source atom mapping')
    coords = {source_key(graph.meta[k]): np.array(tuple(a.pos)) for k, a in graph.atoms.items()}
    metal = np.asarray(repair['selected_site']['xyz_A'], dtype=float)
    direction = coords[source_key(graph.meta[o])] - metal
    direction /= np.linalg.norm(direction)
    origin = coords[source_key(graph.meta[ca])]
    axis = coords[source_key(graph.meta[nb])] - origin
    axis /= np.linalg.norm(axis)
    if not np.isfinite(direction).all() or not np.isfinite(axis).all():
        raise InvalidArtifact('degenerate physical coordinate')
    site = graph.locate(repair['selected_site'])
    mk = graph.key(site.key, repair['selected_site']['atom'])
    if not np.array_equal(coords[source_key(graph.meta[mk])], metal):
        raise InvalidArtifact('metal center does not match pinned source')
    coordinates = {
        'metal': {'units': 'angstrom', 'positive_direction': 'metal_toward_Gln140_carbonyl_O',
                  'axis': direction.tolist(), 'moving_source_keys': ['metal']},
        'peptide': {'units': 'radian', 'positive_direction': 'right_hand_about_CA140_to_CA141',
                    'axis': axis.tolist(), 'origin_A': origin.tolist(),
                    'moving_source_keys': [source_key(graph.meta[k]) for k in moving],
                    'fixed_anchor_keys': [source_key(graph.meta[k]) for k in (ca, nb)]},
    }
    return graph, coords, metal, mk, coordinates


def moved_source(coords, metal, spec, amount):
    moved = {k: v.copy() for k, v in coords.items()}
    newmetal = metal.copy()
    if spec['units'] == 'angstrom':
        newmetal += amount * np.asarray(spec['axis'])
    else:
        origin = np.asarray(spec['origin_A'])
        for key in spec['moving_source_keys']:
            moved[key] = origin + rotation(coords[key] - origin, spec['axis'], amount)
    return moved, newmetal


def physical_velocity(coords, spec):
    if spec['units'] == 'angstrom':
        return {'metal': list(spec['axis'])}
    return {key: np.cross(spec['axis'], coords[key] - spec['origin_A']).tolist()
            for key in spec['moving_source_keys']}


def materialize_motion(repair, coords, metal, velocity):
    """Apply the same graph/cap construction to displaced source coordinates."""
    rows = repair['atom_graph']['source_to_qm']
    if sorted(row['qm_index'] for row in rows) != list(range(1, len(rows) + 1)):
        raise InvalidArtifact('nonbijective QM atom index mapping')
    atoms = [None] * (len(rows) + 1)
    atoms[0] = ('La', *metal.tolist())
    jacobian = np.zeros((len(atoms), 3))
    jacobian[0] = velocity.get('metal', [0, 0, 0])
    cap_indices = []
    for row in rows:
        i = row['qm_index']
        if row['kind'] == 'source':
            key = source_key(row['source'])
            atoms[i] = (row['source']['element'], *coords[key].tolist())
            jacobian[i] = velocity.get(key, [0, 0, 0])
        elif row['kind'] == 'sigma_link_H':
            a, b = [source_key(row[k]) for k in ('retained', 'omitted')]
            x, y = coords[a], coords[b]
            ja, jb = cap_jacobians(x, y, row['length_A'])
            pos = x + row['length_A'] * (y - x) / np.linalg.norm(y - x)
            atoms[i] = ('H', *pos.tolist())
            jacobian[i] = ja @ np.asarray(velocity.get(a, [0, 0, 0])) + jb @ np.asarray(velocity.get(b, [0, 0, 0]))
            cap_indices.append(i)
        else:
            raise InvalidArtifact('unsupported source-to-QM coordinate map')
    if np.any(jacobian[cap_indices] != 0):
        raise InvalidArtifact('approved directions unexpectedly move link-cap anchors')
    return atoms, jacobian


def membership(graph, coords, metal, metal_key):
    """Reapply the established donor selector on the same full physical source."""
    structure = graph.structure.clone()
    for key, original in graph.atoms.items():
        target = structure[0][key[0]][key[1]].find_atom(key[2], original.altloc)
        if not target:
            raise InvalidArtifact('source atom disappeared while mapping displacement')
        pos = metal if key == metal_key else coords[source_key(graph.meta[key])]
        target.pos = gemmi.Position(*pos)
    indexed, _ = generic.index_residues(structure[0])
    original = graph.atoms[metal_key]
    m = structure[0][metal_key[0]][metal_key[1]].find_atom(metal_key[2], original.altloc)
    contacts = generic.collect_contacts(indexed, m, qm_inclusion_cut_A=3.3)
    return {name: sorted((c.chain, c.resnum, c.insertion_code, c.resname, c.atom, c.donor_type)
                         for c in contacts if getattr(c, name))
            for name in ('qualifies', 'included_in_qm')}


def validate_motion(repair, graph, coords, metal, metal_key, spec, amount):
    moved, newmetal = moved_source(coords, metal, spec, amount)
    maxmotion = float(np.linalg.norm(newmetal - metal))
    for key, meta in graph.meta.items():
        if meta['element'] != 'H':
            sk = source_key(meta)
            maxmotion = max(maxmotion, float(np.linalg.norm(moved[sk] - coords[sk])))
    if maxmotion > .05 + 1e-12:
        raise InvalidArtifact('physical heavy-atom displacement exceeds approved 0.05 A trust region')
    if membership(graph, coords, metal, metal_key) != membership(graph, moved, newmetal, metal_key):
        raise InvalidArtifact('displacement changes the frozen typed coordination membership')
    atoms, _ = materialize_motion(repair, moved, newmetal, {})
    return atoms, {'maximum_source_heavy_atom_displacement_A': maxmotion,
                   'coordination_membership': 'unchanged', 'protonation': 'unchanged',
                   'explicit_water_inventory': repair['explicit_water_inventory'],
                   'link_caps_reconstructed': True, 'link_cap_anchor_motion': False}


def validate_representation(name, path):
    repair = read_json(path)
    if repair['protocol_id'] != REPRESENTATIONS[name] or repair.get('explicit_water_inventory') != []:
        raise InvalidArtifact('Stage C only supports the approved dry GGR extended/connected models')
    for metal in ('La', 'Ca'):
        out = repair['outputs'][metal]
        state = input_state(verify(out['input']), gradient=False, tight=False)
        verify(out['xyz'])
        if state['charge'] != out['charge'] or state['multiplicity'] != out['multiplicity']:
            raise InvalidArtifact('representation charge/multiplicity mismatch')
        if out['charge'] != (0 if metal == 'La' else -1):
            raise InvalidArtifact('approved ligand charge -3 not preserved')
    paired(verify(repair['outputs']['La']['xyz']), verify(repair['outputs']['Ca']['xyz']), 0, -1)
    graph, coords, metal, metal_key, specs = physical_model(repair)
    regenerated, _ = materialize_motion(repair, coords, metal, {})
    actual = xyz(verify(repair['outputs']['La']['xyz']))
    if [a[0] for a in regenerated] != [a[0] for a in actual] or not np.allclose(
            [a[1:] for a in regenerated], [a[1:] for a in actual], atol=1e-9, rtol=0):
        raise InvalidArtifact('frozen representation disagrees with physical graph construction')
    return repair, graph, coords, metal, metal_key, specs


def write_task(root, name, repair_path, repair, label, metal, atoms, *, amount=0., coordinate=None, phase='nominal'):
    tid = f'{name}_{label}_{metal}'
    directory = root / tid
    directory.mkdir()
    xp, ip = directory / 'core.xyz', directory / 'endpoint.inp'
    center = label == 'center'
    if center:
        shutil.copyfile(verify(repair['outputs'][metal]['xyz']), xp)
    else:
        lines = [str(len(atoms)), f'{PROTOCOL}; {name}; {label}; {metal}']
        lines += [f'{metal if i == 0 else row[0]} {row[1]:.12f} {row[2]:.12f} {row[3]:.12f}'
                  for i, row in enumerate(atoms)]
        xp.write_text('\n'.join(lines) + '\n')
    charge = repair['outputs'][metal]['charge']
    ip.write_text(f'! {METHOD}{" EnGrad" if center else ""}\n* xyzfile {charge} 1 core.xyz\n')
    settings = {'protocol_id': PROTOCOL, 'representation': record(repair_path),
                'source_structure': repair['source_structure'], 'source_xyz': repair['outputs'][metal]['xyz'],
                'method_input': record(ip), 'coordinates': record(xp), 'metal': metal,
                'charge': charge, 'multiplicity': 1, 'coordinate': coordinate,
                'displacement': amount, 'phase': phase, 'task_type': 'analytic_gradient' if center else 'single_point'}
    return {'task_id': tid, 'input': record(ip), 'xyz': record(xp), 'output_path': str(directory / 'endpoint.out'),
            'engrad_path': str(directory / 'endpoint.engrad') if center else None,
            'metal': metal, 'charge': charge, 'multiplicity': 1, 'representation': name,
            'coordinate': coordinate, 'amplitude': amount, 'phase': phase,
            'task_type': settings['task_type'], 'cache_key': cache_key(settings), 'scientific_settings': settings}


def manifest_header(root, agreement, plan):
    verify(record(agreement)); verify(record(plan))
    implementation = root / 'implementation'
    implementation.mkdir()
    pins = {}
    for name in ('ggr_sensitivity.py', 'ggr_preparation.py', 'affordable_response.py', 'affordable_common.py',
                 'affordable_peptide.py', 'carve_generic.py', 'coordination_policy.py'):
        source = ROOT / 'scripts' / name
        if not source.is_file():
            raise InvalidArtifact(f'missing implementation dependency: {source}')
        destination = implementation / name
        shutil.copyfile(source, destination)
        pins[name] = {'live_source': record(source), 'preserved_copy': record(destination)}
    return {'schema_version': 'alquemia.ggr_sensitivity.v1', 'protocol_id': PROTOCOL,
            'agreement': record(agreement), 'plan': record(plan), 'orca': record(ORCA),
            'implementation': pins, 'execution_policy': {
                'task_runner': record(ROOT / 'scripts/run_orca_task_manifest.py'),
                'runtime_renderer': record(ROOT / 'scripts/render_orca_runtime_input.py')},
            'cost_tracking': {'compute_budget': None, 'wall_time_limit': None},
            'energy_scope': 'isolated_CPCM_endpoint', 'environment_gradient_included': False,
            'conditional_policy': {'nominal_absolute_tolerance_kcal_mol': .02,
                                   'half_absolute_tolerance_kcal_mol': .01,
                                   'relative_tolerance': .05, 'maximum_half_step_tasks': 16,
                                   'normal_tight_bridge_tolerance_kcal_mol': .05}, **UNAVAILABLE}


def prepare(representations, output, agreement, plan, normal_manifest=None):
    root = Path(output).resolve()
    if set(representations) != set(REPRESENTATIONS):
        raise InvalidArtifact('exactly extended and connected representations are required')
    if root.exists():
        raise InvalidArtifact(f'refusing existing workspace: {root}')
    prepared = {name: validate_representation(name, path) for name, path in representations.items()}
    root.mkdir(parents=True)
    result = manifest_header(root, agreement, plan)
    result['normal_manifest'] = record(normal_manifest) if normal_manifest is not None else None
    result.update(tasks=[], representations={}, directions=[], phase='nominal')
    for name, (repair, graph, coords, metal, mk, specs) in prepared.items():
        rp = Path(representations[name]).resolve()
        result['representations'][name] = {'preparation': record(rp), 'normal_outputs': repair.get('tasks', [])}
        for endpoint in ('La', 'Ca'):
            result['tasks'].append(write_task(root, name, rp, repair, 'center', endpoint, None))
        for coordinate, spec in specs.items():
            h = AMPLITUDES[coordinate]
            velocity = physical_velocity(coords, spec)
            _, jacobian = materialize_motion(repair, coords, metal, velocity)
            row = {'representation': name, 'coordinate': coordinate, 'amplitude': h,
                   'physical_coordinate': spec, 'source_jacobian': velocity,
                   'qm_jacobian': jacobian.tolist(), 'units': 'kcal_mol_per_A' if coordinate == 'metal' else 'kcal_mol_per_radian',
                   'moving_cap_energy_validation': False, 'task_ids': [], 'status': 'prepared'}
            try:
                displacements = {sign: validate_motion(repair, graph, coords, metal, mk, spec, sign*h)
                                 for sign in (-1, 1)}
            except InvalidArtifact as exc:
                row.update(status='unsupported', reason=str(exc))
                result['directions'].append(row)
                continue
            for sign, (atoms, validation) in displacements.items():
                for endpoint in ('La', 'Ca'):
                    label = coordinate + ('_plus' if sign == 1 else '_minus')
                    task = write_task(root, name, rp, repair, label, endpoint, atoms,
                                      amount=sign*h, coordinate=coordinate)
                    task['physical_validation'] = validation
                    row['task_ids'].append(task['task_id'])
                    result['tasks'].append(task)
            result['directions'].append(row)
    result['status'] = 'approved_prepared'
    write_new(root / 'manifest.json', result)
    return result


def endpoint_output_validation(task, output, *, require_gradient=False, require_tight=True):
    """Validate the actual Hamiltonian/state/ECP, beyond an energy regex."""
    state = input_state(verify(task['input']), gradient=require_gradient, tight=require_tight)
    if state['charge'] != task['charge'] or state['multiplicity'] != task['multiplicity']:
        raise InvalidArtifact('manifest/input charge or multiplicity differs')
    atoms = xyz(verify(task['xyz']))
    if atoms[0][0] != task['metal'] or any(a[0] not in {'La', 'Ca', 'C', 'H', 'N', 'O'} for a in atoms):
        raise InvalidArtifact('unexpected Stage C element inventory')
    if any(a[0] in {'La', 'Ca'} for a in atoms[1:]):
        raise InvalidArtifact('unexpected additional metal')
    text = Path(output).read_text()
    if not re.search(r'Program Version\s+6\.1\.1\b', text):
        raise InvalidArtifact('actual ORCA version is not approved 6.1.1')
    def integer(pattern, label):
        hits = re.findall(pattern, text, re.I)
        if len(hits) != 1:
            raise InvalidArtifact(f'output {label} absent/ambiguous')
        return int(hits[0])
    charge = integer(r'Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)', 'charge')
    multiplicity = integer(r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)', 'multiplicity')
    electrons = integer(r'Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)', 'electron count')
    ecp = re.findall(r'Type\s+(\w+)\s+ECP\s+(\S+)\s+\(replacing\s+(\d+)\s+core electrons', text)
    expected_ecp = [('La', 'Def2-ECP', '46')] if task['metal'] == 'La' else []
    if ecp != expected_ecp:
        raise InvalidArtifact('unexpected native ECP convention')
    expected_electrons = sum(gemmi.Element(a[0]).atomic_number for a in atoms) - task['charge'] - (46 if task['metal'] == 'La' else 0)
    if (charge, multiplicity, electrons) != (task['charge'], 1, expected_electrons) or electrons % 2:
        raise InvalidArtifact('output charge/multiplicity/ECP electron parity mismatch')
    for pattern, label in [(r'Surface type\s+\.{2,}\s+GAUSSIAN VDW', 'Gaussian VDW CPCM'),
                           (r'CPCM SOLVATION MODEL', 'CPCM'), (r'DFTD4', 'native D4'),
                           (r'gCP correction\s+[-+0-9.]', 'native gCP')]:
        if not re.search(pattern, text, re.I):
            raise InvalidArtifact(f'actual output missing {label}')
    if re.search(r'numerical (?:gradient|differentiation)', text, re.I):
        raise InvalidArtifact('numerical DFT differentiation is unsupported')
    gradient_components = None
    if require_gradient:
        patterns = {
            'analytic_scf': r'ORCA SCF GRADIENT CALCULATION',
            'cpcm': r'CPCM gradient\s+\.{2,}\s+done',
            'native_dispersion': r'DISPERSION GRADIENT',
            'native_gcp': r'gCP correction\s+\.{2,}\s+done',
            'total_cartesian': r'CARTESIAN GRADIENT',
        }
        if task['metal'] == 'La':
            patterns['native_ecp'] = r'ECP gradient\s+\(SHARK\)\s+\.{2,}\s+done'
        gradient_components = {name: bool(re.search(pattern, text, re.I)) for name, pattern in patterns.items()}
        if not all(gradient_components.values()):
            raise InvalidArtifact('actual output lacks a required native analytic-gradient component')
    return {'orca_version': '6.1.1', 'charge': charge, 'multiplicity': multiplicity,
            'explicit_electrons': electrons, 'ecp_core_electrons': 46 if task['metal'] == 'La' else 0,
            'total_energy_includes_native_D4_gCP_CPCM': True,
            'analytic_gradient_components_reported': gradient_components,
            'composite_gradient_status': 'requires_directional_energy_consistency' if require_gradient else 'not_a_gradient'}


def executed(manifest_path):
    """Read completed tasks only after the existing runner verifies full receipts."""
    from run_orca_task_manifest import load_manifest_tasks, _completed_attempt_is_valid
    m, tasks = load_manifest_tasks(Path(manifest_path))
    verify(m['agreement']); verify(m['orca'])
    if 'plan' in m:
        verify(m['plan'])
    for dependency in m.get('implementation', {}).values():
        if isinstance(dependency, dict) and 'preserved_copy' in dependency:
            verify(dependency['preserved_copy'])
    for pin in m['execution_policy'].values():
        verify(pin)
    result = {}
    for task in tasks:
        key = task['task_id']
        verify({'path': str(task['input']), 'sha256': task['input_sha256']})
        verify({'path': str(task['xyz']), 'sha256': task['xyz_sha256']})
        output = task['output']
        receipt = Path(str(output) + '.execution.json')
        if not output.exists() or not receipt.exists():
            result[key] = {'status': 'missing', 'reason': 'output or execution receipt unavailable'}
            continue
        if not _completed_attempt_is_valid(receipt, output, manifest_sha256=digest(manifest_path), task=task,
                runner_identity=m['execution_policy']['task_runner'], runtime_renderer_identity=m['execution_policy']['runtime_renderer']):
            result[key] = {'status': 'invalid', 'reason': 'execution receipt failed provenance/completion validation',
                           'receipt': record(receipt), 'output': record(output)}
            continue
        try:
            receipt_data = read_json(receipt)
            if receipt_data.get('orca_executable') != m['orca']:
                raise InvalidArtifact('execution used a different ORCA executable')
            for artifact in receipt_data['artifacts'].values():
                verify(artifact)
            result[key] = {'status': 'complete', 'energy_hartree': energy(output),
                           'receipt': record(receipt), 'output': record(output)}
        except InvalidArtifact as exc:
            result[key] = {'status': 'invalid', 'reason': str(exc)}
    return m, result


def consistency(eplus, eminus, ecenter, derivative, h, *, half=False):
    if h <= 0 or not all(math.isfinite(x) for x in (eplus, eminus, ecenter, derivative, h)):
        raise InvalidArtifact('nonfinite directional energy/gradient or invalid amplitude')
    odd = (eplus - eminus) * HA_TO_KCAL / 2
    even = ((eplus - ecenter) + (eminus - ecenter)) * HA_TO_KCAL / 2
    prediction = h * derivative
    tolerance = max(.01 if half else .02, .05 * abs(prediction))
    residual = odd - prediction
    return {'odd_energy_kcal_mol': odd, 'predicted_odd_energy_kcal_mol': prediction,
            'residual_kcal_mol': residual, 'residual_per_coordinate_unit': residual/h,
            'even_curvature_sensitive_energy_kcal_mol': even, 'tolerance_kcal_mol': tolerance,
            'projected_gradient': derivative, 'amplitude': h, 'pass': abs(residual) <= tolerance}


def normal_bridge(m, energies, normal_manifest):
    if normal_manifest is None:
        return {'status': 'unavailable', 'reason': 'provide --normal-manifest with executed Stage A task receipts'}
    nm, ne = executed(normal_manifest)
    result = {'normal_manifest': record(normal_manifest), 'representations': {}}
    for name, rep in m['representations'].items():
        repair = read_json(verify(rep['preparation']))
        changes = {}
        for metal in ('La', 'Ca'):
            wanted = repair['outputs'][metal]['xyz']['sha256']
            matches = [t for t in nm['tasks'] if t['xyz']['sha256'] == wanted]
            if len(matches) != 1 or ne[matches[0]['task_id']]['status'] != 'complete':
                changes[metal] = {'status': 'unavailable'}
                continue
            normal = matches[0]
            state = input_state(verify(normal['input']), gradient=False, tight=False)
            if state['charge'] != repair['outputs'][metal]['charge']:
                raise InvalidArtifact('normal-SCF bridge endpoint charge differs')
            endpoint_output_validation({**normal, 'metal': metal, 'charge': state['charge'], 'multiplicity': 1},
                                       verify(ne[normal['task_id']]['output']), require_tight=False)
            tight = energies.get(f'{name}_center_{metal}', {})
            if tight.get('status') != 'complete':
                changes[metal] = {'status': 'unavailable'}
                continue
            delta = (tight['energy_hartree'] - ne[normal['task_id']]['energy_hartree']) * HA_TO_KCAL
            changes[metal] = {'status': 'complete', 'tight_minus_normal_kcal_mol': delta,
                              'normal_output': ne[normal['task_id']]['output'], 'pass': abs(delta) <= .05}
        if all(changes[x]['status'] == 'complete' for x in ('La', 'Ca')):
            dr = changes['Ca']['tight_minus_normal_kcal_mol'] - changes['La']['tight_minus_normal_kcal_mol']
            changes['R'] = {'tight_minus_normal_kcal_mol': dr, 'pass': abs(dr) <= .05}
            changes['interchangeable_at_declared_scale'] = all(changes[x]['pass'] for x in ('La', 'Ca', 'R'))
        result['representations'][name] = changes
    return result


def collect(manifest_path, normal_manifest=None):
    m, energies = executed(manifest_path)
    pinned_normal = m.get('normal_manifest')
    if pinned_normal is not None:
        pinned_path = verify(pinned_normal)
        if normal_manifest is not None and record(normal_manifest) != pinned_normal:
            raise InvalidArtifact('normal-SCF override differs from the frozen Stage A manifest')
        normal_manifest = pinned_path
    half = m['phase'] == 'half'
    center_manifest = verify(m['nominal_manifest']) if half else Path(manifest_path)
    cm, centers = executed(center_manifest) if half else (m, energies)
    gradients = {}
    for task in cm['tasks']:
        if task['task_type'] != 'analytic_gradient':
            continue
        key = task['task_id']
        if centers[key]['status'] != 'complete':
            continue
        try:
            validation = endpoint_output_validation(task, verify(centers[key]['output']), require_gradient=True)
            rp = verify(cm['representations'][task['representation']]['preparation'])
            gradients[key] = extract(task['engrad_path'], verify(centers[key]['output']),
                                     verify(task['input']), verify(task['xyz']), rp)
            gradients[key]['state_validation'] = validation
        except (InvalidArtifact, OSError, ValueError) as exc:
            centers[key] = {**centers[key], 'status': 'invalid', 'reason': str(exc)}
    for task in m['tasks']:
        key = task['task_id']
        if task['task_type'] == 'analytic_gradient' or energies[key]['status'] != 'complete':
            continue
        try:
            energies[key]['state_validation'] = endpoint_output_validation(task, verify(energies[key]['output']))
        except (InvalidArtifact, OSError, ValueError) as exc:
            energies[key] = {**energies[key], 'status': 'invalid', 'reason': str(exc)}
    comparisons = []
    taskmap = {t['task_id']: t for t in m['tasks']}
    for direction in m['directions']:
        row = {'representation': direction['representation'], 'coordinate': direction['coordinate'],
               'phase': m['phase'], 'amplitude': direction['amplitude'], 'units': direction['units'],
               'status': direction['status'], 'half_step_eligible': False, **UNAVAILABLE}
        if direction['status'] != 'prepared':
            row['reason'] = direction.get('reason')
            comparisons.append(row)
            continue
        name = direction['representation']
        keys = {metal: f'{name}_center_{metal}' for metal in ('La', 'Ca')}
        if any(keys[x] not in gradients for x in keys) or any(energies[t]['status'] != 'complete' for t in direction['task_ids']):
            row.update(status='unavailable', reason='missing/invalid center gradient or displaced endpoint')
            comparisons.append(row)
            continue
        derivatives, endpoint_checks, values = {}, {}, {}
        for metal in ('La', 'Ca'):
            gradient = gradients[keys[metal]]
            mapped = gradient['mapped_source_gradient_kcal_mol_per_A']
            velocity = direction['source_jacobian']
            if set(velocity) - set(mapped):
                raise InvalidArtifact('physical displacement includes an unmapped source atom')
            derivatives[metal] = float(sum(np.dot(mapped[key], value) for key, value in velocity.items()))
            qm_projection = float(np.sum(np.asarray(gradient['gradient_kcal_mol_per_A']) * direction['qm_jacobian']))
            if not math.isclose(derivatives[metal], qm_projection, abs_tol=1e-7, rel_tol=1e-10):
                raise InvalidArtifact('source and QM Jacobian projections differ')
            subset = [taskmap[t] for t in direction['task_ids'] if taskmap[t]['metal'] == metal]
            if len(subset) != 2:
                raise InvalidArtifact('direction endpoint multiplicity differs from approved pair')
            plus = next(t for t in subset if t['amplitude'] > 0)
            minus = next(t for t in subset if t['amplitude'] < 0)
            values[metal] = (energies[plus['task_id']]['energy_hartree'], energies[minus['task_id']]['energy_hartree'],
                             centers[keys[metal]]['energy_hartree'])
            endpoint_checks[metal] = consistency(*values[metal], derivatives[metal], direction['amplitude'], half=half)
        rvalues = [ca - la for ca, la in zip(values['Ca'], values['La'])]
        endpoint_checks['R'] = consistency(*rvalues, derivatives['Ca'] - derivatives['La'], direction['amplitude'], half=half)
        passed = all(v['pass'] for v in endpoint_checks.values())
        paired_map = {key: (np.asarray(gradients[keys['Ca']]['mapped_source_gradient_kcal_mol_per_A'][key]) - value).tolist()
                      for key, value in gradients[keys['La']]['mapped_source_gradient_kcal_mol_per_A'].items()}
        row.update(status='pass' if passed else 'failed_consistency', checks=endpoint_checks,
                   half_step_eligible=not half and not passed,
                   physical_source_grad_E_Ca_minus_E_La_kcal_mol_per_A=paired_map,
                   interpretation='directional_energy_consistency_only; no moving-cap energy validation')
        comparisons.append(row)
    return {'schema_version': 'alquemia.ggr_sensitivity_comparison.v1', 'manifest': record(manifest_path),
            'collection_implementation': {'collector': record(__file__),
                                          'response_adapter': record(Path(extract.__code__.co_filename))},
            'phase': m['phase'], 'energies': energies, 'gradients': gradients, 'comparisons': comparisons,
            'normal_tight_bridge': normal_bridge(cm, centers, normal_manifest),
            'half_step_blocks': [{'representation': r['representation'], 'coordinate': r['coordinate']}
                                 for r in comparisons if r['half_step_eligible']],
            'numerical_cause': 'unresolved; half-step pass alone does not identify the cause',
            'absolute_classification': 'uncalibrated_protocol', **UNAVAILABLE}


def prepare_half(manifest_path, comparison_path, output):
    m = read_json(manifest_path)
    supplied = read_json(comparison_path)
    if m.get('phase') != 'nominal' or supplied.get('manifest') != record(manifest_path):
        raise InvalidArtifact('half-step requires the original nominal manifest and its comparison')
    # Recompute selection from pinned real outputs: an edited status cannot widen scope.
    comparison = collect(manifest_path)
    if supplied['half_step_blocks'] != comparison['half_step_blocks']:
        raise InvalidArtifact('half-step selection disagrees with current verified endpoint evidence')
    selected = {(x['representation'], x['coordinate']) for x in comparison['half_step_blocks']}
    if not selected:
        return {'status': 'no_failed_valid_blocks', 'tasks': []}
    root = Path(output).resolve()
    if root.exists():
        raise InvalidArtifact(f'refusing existing half-step workspace: {root}')
    root.mkdir(parents=True)
    result = manifest_header(root, verify(m['agreement']), verify(m['plan']))
    result.update(phase='half', nominal_manifest=record(manifest_path), trigger_comparison=record(comparison_path),
                  representations=m['representations'], normal_manifest=m.get('normal_manifest'),
                  tasks=[], directions=[], status='approved_conditional_prepared')
    for original in m['directions']:
        name, coordinate = original['representation'], original['coordinate']
        if (name, coordinate) not in selected:
            continue
        rp = verify(m['representations'][name]['preparation'])
        repair, graph, coords, metal, mk, specs = validate_representation(name, rp)
        row = copy.deepcopy(original)
        row['amplitude'] /= 2
        row['task_ids'] = []
        h = row['amplitude']
        displacements = {sign: validate_motion(repair, graph, coords, metal, mk, specs[coordinate], sign*h)
                         for sign in (-1, 1)}
        for sign, (atoms, validation) in displacements.items():
            for endpoint in ('La', 'Ca'):
                label = coordinate + ('_half_plus' if sign == 1 else '_half_minus')
                task = write_task(root, name, rp, repair, label, endpoint, atoms,
                                  amount=sign*h, coordinate=coordinate, phase='half')
                task['physical_validation'] = validation
                row['task_ids'].append(task['task_id'])
                result['tasks'].append(task)
        result['directions'].append(row)
    if len(result['tasks']) > 16:
        raise InvalidArtifact('conditional task list exceeds approved named blocks')
    write_new(root / 'manifest.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='operation', required=True)
    p = commands.add_parser('prepare')
    p.add_argument('--representation', action='append', required=True, metavar='NAME=MANIFEST')
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--agreement', type=Path, required=True)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--normal-manifest', type=Path, required=True)
    p = commands.add_parser('collect')
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--normal-manifest', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p = commands.add_parser('prepare-half')
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--comparison', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.operation == 'prepare':
        entries = [value.split('=', 1) for value in args.representation]
        if any(len(e) != 2 for e in entries) or len({e[0] for e in entries}) != len(entries):
            parser.error('representations must be unique NAME=MANIFEST entries')
        result = prepare(dict(entries), args.output_dir, args.agreement, args.plan, args.normal_manifest)
    elif args.operation == 'collect':
        result = collect(args.manifest, args.normal_manifest)
        write_new(args.output, result)
    else:
        result = prepare_half(args.manifest, args.comparison, args.output_dir)
    print(f'{args.operation}: {len(result.get("tasks", []))} prepared tasks; {result.get("status", "collected")}; response_model_not_validated')


if __name__ == '__main__':
    main()
