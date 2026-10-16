#!/usr/bin/env python3
"""Fixed LBT3 sources plus one rigid actual H2O; geometry only, Slurm required."""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time

PROTOCOL = 'lbt3_fixed_source_one_water_geometry_v1'
PLAN_SHA = 'f44bfdcd0bcee76dcf4336c59efe29a40b3314f2c8d35d230534b691b67eb47d'
PREPARATION_SHA = '649a12d4c78a72058c7446184a67cf252d3b37536c76beb700685713d3e881ec'
WATER_SHA = '523ad04fd4c3c33162302cf50296b08b423ac46a0e8d89fb1a04e1016ca49523'
SOURCES = ('7CCO', '7CCN')
METALS = (('La', 1), ('Dy', 6))
POSES = (('inner_plane0', 2.6, 0), ('inner_plane90', 2.6, 90),
         ('outer_plane0', 4.5, 0), ('outer_plane90', 4.5, 90))
SETTINGS = {'polar_angles_degree': [60, 75], 'azimuths_degree': list(range(0, 360, 30)),
            'inner_O_radius_A': 2.6, 'outer_O_radius_A': 4.5,
            'water_plane_angles_degree': [0, 90], 'frame_minimum_norm': 1e-10,
            'roundoff_A': 1e-10, 'water_shape_tolerance_A': 1e-10,
            'cutoffs_A': {'O_peptide_heavy': 2.0, 'O_peptide_H': 1.2,
                          'H_peptide_heavy': 1.2, 'H_peptide_H': .8, 'H_metal': 2.0},
            'selection': 'max minimum(distance/cutoff) across four poses; first exact tie',
            'original_atoms_fixed': 256, 'water_atoms_added': 3, 'energy_calls': 0}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(value):
    path = Path(value['path']).resolve(strict=True)
    require(path.is_file() and pin(path)['sha256'] == value['sha256'],
            'missing/changed pinned input: ' + str(value['path']))
    return path


def read(value):
    return json.loads(verify(value).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    rows = [line.split() for line in lines[2:] if line.strip()]
    require(len(rows) == int(lines[0]) and all(len(r) == 4 for r in rows), 'invalid XYZ inventory')
    return [(r[0], *[float(v) for v in r[1:]]) for r in rows]


def write_xyz(path, rows):
    with Path(path).open('x') as handle:
        handle.write(str(len(rows)) + '\nModeled one-water placement; original source fixed\n')
        for z, *position in rows:
            handle.write(z + ' ' + ' '.join(format(v, '.17g') for v in position) + '\n')


def inputs(config):
    cfg = json.loads(Path(config).read_text())
    require(cfg['implementation']['sha256'] == pin(__file__)['sha256'], 'implementation hash differs')
    verify(cfg['implementation'])
    for key, expected in (('plan', PLAN_SHA), ('preparation_result', PREPARATION_SHA),
                          ('water_template', WATER_SHA)):
        require(cfg[key]['sha256'] == expected, key + ' differs from declared exact source')
        verify(cfg[key])
    for value in cfg.get('pins', []):
        verify(value)
    repo = Path(cfg['repository']).resolve(strict=True)
    require(set(cfg['dependencies']) == {'hydration_basin_coordinates.py', 'affordable_common.py'},
            'exact two geometry dependencies required')
    for name, value in cfg['dependencies'].items():
        verify(value)
        require(pin(repo / 'scripts' / name)['sha256'] == value['sha256'], 'import dependency differs')
    sys.path.insert(0, str(repo / 'scripts'))
    import numpy as np
    import hydration_basin_coordinates as geometry
    require(pin(geometry.__file__)['sha256'] == cfg['dependencies']['hydration_basin_coordinates.py']['sha256'],
            'loaded rigid-water helper differs')
    import affordable_common
    require(pin(affordable_common.__file__)['sha256'] == cfg['dependencies']['affordable_common.py']['sha256'],
            'loaded helper dependency differs')
    prior = read(cfg['preparation_result'])
    require(prior['status'] == 'complete' and prior['protocol_id'] ==
            'lbt3_model1_whole17_pH6_ff19sb_terminal_completion_LaDy_v1', 'original preparation unavailable')
    for key in ('config', 'plan', 'implementation', 'legacy_builder'):
        verify(prior[key])
    for value in prior['dependencies'].values():
        verify(value)
    require([r['source_id'] for r in prior['rows']] == list(SOURCES), 'source scope/order differs')
    template = read_xyz(verify(cfg['water_template']))
    require([r[0] for r in template] == ['O', 'H', 'H'], 'complete pinned neutral H2O required')
    require(np.isfinite(np.asarray([r[1:] for r in template])).all(), 'nonfinite H2O geometry')
    prepared = []
    for row in prior['rows']:
        p = read(row['preparation'])
        require(row['status'] == 'prepared' and p['status'] == 'prepared' and
                p['source_id'] == row['source_id'] and p['tasks'] == row['tasks'], 'source receipt differs')
        require(p['sequence'] == 'FIDTNNDGWIEGDELLA' and p['model'] == 1 and
                p['atoms'] == 256 and p['protein_atoms'] == 255 and p['charge'] == -2 and
                p['protein_charge'] == -5 and p['waters'] == p['spectators'] == p['caps'] == [],
                'source chemistry/inventory differs')
        require(p['deposited_H_coordinates_preserved'] and p['deposited_heavy_coordinates_preserved'],
                'deposited atom preservation missing')
        for value in p['source_pins'].values():
            verify(value)
        for key in ('native_peptide_system', 'forcefield', 'hydrogen_audit', 'completed_peptide'):
            verify(p[key])
        atoms = read(p['physical_atoms']); native = read(p['native_peptide_atoms'])
        bonds = read(p['covalent_bonds']); ids = [a['id'] for a in atoms]
        require(len(atoms) == len(set(ids)) == 256 and ids[-1] == 'metal' and
                atoms[-1]['element'] == 'M', 'physical particle map differs')
        require(len(native) == 255 and all(a['id'] == b['id'] and a['element'] == b['element'] and
                a['xyz_A'] == b['xyz_A'] for a, b in zip(atoms[:-1], native)), 'native peptide map differs')
        require(all(a['element'] in ('C', 'N', 'O', 'H') for a in atoms[:-1]), 'unexpected source species')
        x = np.asarray([a['xyz_A'] for a in atoms], dtype=float)
        require(x.shape == (256, 3) and np.isfinite(x).all(), 'invalid source coordinates')
        for b in bonds:
            i, j = b['native_particle_indices']
            require(0 <= i < 255 and 0 <= j < 255 and ids[i] == b['atom_a_id'] and
                    ids[j] == b['atom_b_id'], 'native bond/source mapping differs')
        for metal, mult in METALS:
            task = p['endpoints']['whole'][metal]
            require(task in row['tasks'] and task['source_id'] == row['source_id'] and
                    task['metal'] == metal and task['charge'] == -2 and task['multiplicity'] == mult and
                    task['spin_multiplicity'] == mult and task['metal_index'] == 255 and task['atom_count'] == 256
                    and task['physical_source'] == p['physical_atoms'], 'paired source state differs')
            expected = [(metal if a['element'] == 'M' else a['element'], *a['xyz_A']) for a in atoms]
            require(read_xyz(verify(task['xyz'])) == expected, 'paired XYZ not exact original physical source')
        for name in ('A/5//OD1', 'A/5//CG'):
            require(name in ids and atoms[ids.index(name)]['resname'] == 'ASN', 'Asn5 frame mapping missing')
        prepared.append({'receipt': row, 'preparation': p, 'atoms': atoms, 'bonds': bonds})
    return cfg, prepared, template, geometry.WaterCoordinates


def unit(v):
    import numpy as np
    norm = float(np.linalg.norm(v))
    require(norm > SETTINGS['frame_minimum_norm'], 'degenerate declared source/template frame')
    return v / norm


def contacts(atoms, position):
    import numpy as np
    original = np.asarray([a['xyz_A'] for a in atoms])
    cutoffs = SETTINGS['cutoffs_A']; values = []; violations = []; nearest = []
    for k, element in enumerate(('O', 'H', 'H')):
        local = []
        for i, a in enumerate(atoms):
            distance = float(np.linalg.norm(position[256 + k] - original[i]))
            if i == 255 and element == 'O':
                continue  # Exact radial check is separate; no ray score contribution.
            key = 'H_metal' if i == 255 else element + ('_peptide_H' if a['element'] == 'H' else '_peptide_heavy')
            cutoff = cutoffs[key]
            row = {'water_index': 256 + k, 'original_index': i, 'original_id': a['id'],
                   'distance_A': distance, 'cutoff_A': cutoff, 'pair_class': key,
                   'normalized_clearance': distance / cutoff}
            values.append(distance / cutoff); local.append(row)
            if distance < cutoff - SETTINGS['roundoff_A']:
                violations.append(row)
        nearest.extend(sorted(local, key=lambda r: (r['distance_A'], r['original_index']))[:4])
    return {'pass': not violations, 'minimum_normalized_clearance': min(values),
            'nearest_contacts_by_water_atom': nearest, 'violations': violations,
            'inherited_source_contacts_unchanged': True}


def prepare_source(source, template, template_pin, out, WaterCoordinates):
    import numpy as np
    from scipy.spatial.transform import Rotation
    p = source['preparation']; atoms = source['atoms']; ids = [a['id'] for a in atoms]
    x = np.asarray([a['xyz_A'] for a in atoms]); t = np.asarray([a[1:] for a in template])
    metal = x[255]; od = x[ids.index('A/5//OD1')]; cg = x[ids.index('A/5//CG')]
    e0 = unit(od - metal); transverse = cg - od
    e1 = unit(transverse - np.dot(transverse, e0) * e0); e2 = np.cross(e0, e1)
    # The actual template need not be exactly symmetric: retain all components.
    h1, h2 = t[1] - t[0], t[2] - t[0]
    tb = unit(h1 + h2); diff = h1 - h2
    tt = unit(diff - np.dot(diff, tb) * tb); tn = np.cross(tb, tt)
    template_frame = np.column_stack((tt, tn, tb))
    require(np.linalg.det(template_frame) > 0, 'improper template frame')
    original_rows = [(a['element'], *a['xyz_A']) for a in atoms]
    group = {'role': 'variable', 'indices': [256, 257, 258],
             'oxygen_index': 256, 'hydrogen_indices': [257, 258]}
    helper = WaterCoordinates(original_rows + template, [group])
    require(helper.fixed == list(range(256)) and helper.mobile == [256, 257, 258], 'rigid helper scope differs')
    shape0 = np.linalg.norm(t[:, None, :] - t[None, :, :], axis=2)
    candidates = []; best = None
    for theta in SETTINGS['polar_angles_degree']:
        for phi in SETTINGS['azimuths_degree']:
            th, ph = np.deg2rad([theta, phi])
            ray = unit(np.cos(th) * e0 + np.sin(th) * (np.cos(ph) * e1 + np.sin(ph) * e2))
            plane0 = unit(e0 - np.dot(e0, ray) * ray)
            candidate = {'candidate_index': len(candidates), 'theta_degree': theta, 'phi_degree': phi,
                         'ray_unit': ray.tolist(), 'poses': []}
            geometries = []
            for name, radius, angle in POSES:
                tangent = plane0 if angle == 0 else np.cross(ray, plane0)
                frame = np.column_stack((tangent, np.cross(ray, tangent), ray))
                rot = frame @ template_frame.T
                require(np.allclose(rot.T @ rot, np.eye(3), rtol=0, atol=1e-12) and
                        abs(float(np.linalg.det(rot)) - 1) <= 1e-12, 'water rotation is not proper')
                center = helper.centers[0]
                oxygen = metal + radius * ray
                translation = oxygen - (rot @ (t[0] - center) + center)
                q = np.r_[translation, Rotation.from_matrix(rot).as_rotvec()]
                y = helper.positions(q)
                require(np.array_equal(y[:256], x), 'original source atom moved')
                shape = np.linalg.norm(y[256:, None, :] - y[None, 256:, :], axis=2)
                error = float(np.max(np.abs(shape - shape0)))
                require(error <= SETTINGS['water_shape_tolerance_A'], 'actual water internal geometry changed')
                radial_error = abs(float(np.linalg.norm(y[256] - metal)) - radius)
                require(radial_error <= SETTINGS['roundoff_A'], 'declared oxygen radius differs')
                check = contacts(atoms, y)
                data = {'pose_id': name, 'shell': 'inner' if radius == 2.6 else 'outer',
                        'O_radius_A': radius, 'plane_angle_degree': angle, 'admitted': check['pass'],
                        'contact_check': check, 'water_shape_max_error_A': error,
                        'oxygen_radial_error_A': radial_error, 'original256_exactly_unchanged': True,
                        'rigid_coordinate_q': q.tolist(), 'row_vector_rotation': rot.T.tolist(),
                        'water_COM_translation_A': translation.tolist()}
                candidate['poses'].append(data); geometries.append(y)
            candidate['minimum_normalized_clearance'] = min(
                a['contact_check']['minimum_normalized_clearance'] for a in candidate['poses'])
            candidate['all_four_admitted'] = all(a['admitted'] for a in candidate['poses'])
            candidates.append(candidate)
            if best is None or candidate['minimum_normalized_clearance'] > best[0]['minimum_normalized_clearance']:
                best = (candidate, geometries)
    require(len(candidates) == 24, 'bounded candidate catalogue differs')
    out.mkdir(parents=True, exist_ok=False)
    chosen, geometries = best
    frame_record = {'source_id': p['source_id'], 'source_preparation': source['receipt']['preparation'],
                    'original_physical_source': p['physical_atoms'],
                    'frame_atom_ids': ['metal', 'A/5//OD1', 'A/5//CG'],
                    'frame_physical_indices': [255, ids.index('A/5//OD1'), ids.index('A/5//CG')],
                    'e0': e0.tolist(), 'e1': e1.tolist(), 'e2': e2.tolist(),
                    'selected_candidate_index': chosen['candidate_index'],
                    'candidate_order': 'theta60 then75; azimuth ascending0..330 by30',
                    'selection_uses_only_geometry': True, 'candidates': candidates}
    write(out / 'PLACEMENT.json', frame_record)
    mapping = [{'physical_index': i, 'physical_id': a['id'], 'kind': 'original_source',
                'source_index': i, 'source_atom': a} for i, a in enumerate(atoms)]
    for k, element in enumerate(('O', 'H', 'H')):
        mapping.append({'physical_index': 256 + k, 'physical_id': ('W/1//O', 'W/1//H1', 'W/1//H2')[k],
                        'kind': 'modeled_water', 'element': element, 'template_index': k,
                        'template': template_pin, 'deposited_LBT3_water': False})
    graph = {'original_peptide_bonds': p['covalent_bonds'], 'original_peptide_system': p['native_peptide_system'],
             'water_bonds': [[256, 257], [256, 258]], 'metal_covalent_bonds': [],
             'whole_peptide_bonds': source['bonds'], 'water_internal_pair_distances_A': shape0.tolist(),
             'water_formal_charge': 0, 'water_template': template_pin}
    write(out / 'SOURCE_MAP.json', mapping); write(out / 'GRAPH.json', graph)
    tasks = []; pose_rows = []; atomic_numbers = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'La': 57, 'Dy': 66}
    for data, y in zip(chosen['poses'], geometries):
        directory = out / data['pose_id']; directory.mkdir()
        physical = [dict(a) for a in atoms] + [
            {'id': mapping[i]['physical_id'], 'element': element, 'xyz_A': y[i].tolist(),
             'kind': 'modeled_water', 'template_index': i - 256, 'water_formal_charge': 0}
            for i, element in zip((256, 257, 258), ('O', 'H', 'H'))]
        write(directory / 'PHYSICAL.json', physical)
        pose = dict(data, source_id=p['source_id'], selected_candidate_index=chosen['candidate_index'],
                    physical_source=pin(directory / 'PHYSICAL.json'), source_map=pin(out / 'SOURCE_MAP.json'),
                    graph=pin(out / 'GRAPH.json'), placement=pin(out / 'PLACEMENT.json'))
        for element, mult in METALS:
            symbols = [element if a['element'] == 'M' else a['element'] for a in physical]
            rows = [(z, *v) for z, v in zip(symbols, y)]
            xp = directory / (element + '.xyz'); write_xyz(xp, rows)
            require(read_xyz(xp) == rows, 'XYZ serialization changed original or water coordinates')
            electrons = sum(atomic_numbers[z] for z in symbols) + 2
            require((electrons - (mult - 1)) % 2 == 0, 'paired electron/spin parity differs')
            tasks.append({'task_id': p['source_id'] + '_' + data['pose_id'] + '_' + element,
                          'source_id': p['source_id'], 'pose_id': data['pose_id'], 'shell': data['shell'],
                          'plane_angle_degree': data['plane_angle_degree'], 'representation': 'whole_plus_one_water',
                          'admitted': data['admitted'], 'status': 'prepared' if data['admitted'] else 'geometry_rejected',
                          'xyz': pin(xp), 'metal': element, 'metal_index': 255, 'charge': -2,
                          'multiplicity': mult, 'spin_multiplicity': mult, 'physical_multiplicity': mult,
                          'physical_spin_2S': mult - 1, 'atom_count': 259, 'all_electron_count': electrons,
                          'water_count': 1, 'water_indices': [256, 257, 258], 'oxidation_state': 3,
                          'physical_source': pose['physical_source'], 'source_map': pose['source_map'],
                          'graph': pose['graph'], 'placement': pose['placement'],
                          'original_source': p['physical_atoms'], 'original_endpoint': p['endpoints']['whole'][element],
                          'water_template': template_pin})
        write(directory / 'POSE.json', pose); pose_rows.append(dict(pose, result=pin(directory / 'POSE.json')))
    result = {'source_id': p['source_id'], 'status': 'prepared', 'source_preparation': source['receipt']['preparation'],
              'original_physical_source': p['physical_atoms'], 'source_map': pin(out / 'SOURCE_MAP.json'),
              'graph': pin(out / 'GRAPH.json'), 'placement': pin(out / 'PLACEMENT.json'),
              'selected_candidate_index': chosen['candidate_index'], 'poses': pose_rows, 'tasks': tasks,
              'admitted_poses': sum(r['admitted'] for r in pose_rows), 'expected_poses': 4,
              'all_original256_exact': True, 'paired_coordinates_identical_except_target_symbol': True,
              'terminal_chemistry': p['terminal_chemistry'], 'water_sites_are_modeled': True}
    write(out / 'RESULT.json', result)
    return dict(result, result=pin(out / 'RESULT.json'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(bool(os.environ.get('SLURM_JOB_ID')), 'geometry calculations require Slurm')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    result = {'protocol_id': PROTOCOL, 'status': 'unavailable', 'config': pin(args.config),
              'implementation': pin(__file__), 'settings': SETTINGS, 'rows': [], 'tasks': [],
              'MACE_calls': 0, 'GFN2_calls': 0, 'quantum_calls': 0, 'minimization_calls': 0,
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'expected_sources': list(SOURCES),
              'expected_pose_order': [p[0] for p in POSES], 'expected_task_count': 16}
    try:
        cfg, prepared, template, WaterCoordinates = inputs(args.config)
        result.update(plan=cfg['plan'], preparation_result=cfg['preparation_result'],
                      water_template=cfg['water_template'], dependencies=cfg['dependencies'])
        for source in prepared:
            row = prepare_source(source, template, cfg['water_template'], out / source['preparation']['source_id'],
                                 WaterCoordinates)
            result['rows'].append(row); result['tasks'].extend(row['tasks'])
            print(json.dumps({'source_id': row['source_id'], 'admitted_poses': row['admitted_poses']}), flush=True)
        require(len(result['tasks']) == 16 and len({r['task_id'] for r in result['tasks']}) == 16,
                'complete sixteen paired slots required')
        admitted = sum(t['admitted'] for t in result['tasks'])
        result.update(status='complete', admitted_tasks=admitted, geometry_rejected_tasks=16 - admitted,
                      all_tasks_admitted=admitted == 16, placement_candidates_evaluated=48,
                      physical_poses_evaluated=192, model_evaluations=0,
                      raw_minimum_across_water_inventories_permitted=False,
                      outer_shell_is_bulk_reference=False,
                      occupancy_population_or_affinity_available=False)
    except Exception as exc:
        result.update(status='execution_failure', reason=str(exc), exception=type(exc).__name__)
    finally:
        result.update(wall_seconds=time.monotonic() - start,
                      software={p: importlib.metadata.version(p) for p in ('numpy', 'scipy')})
        write(out / 'RESULT.json', result)
    if result['status'] != 'complete':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
