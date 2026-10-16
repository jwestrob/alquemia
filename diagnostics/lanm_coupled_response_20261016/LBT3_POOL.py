#!/usr/bin/env python3
"""Frozen LBT3 terminal common pool; energy-only cross evaluation, Slurm only.

CLI: LBT3_POOL.py {run,report} --config CONFIG --output FRESH
Config: repository, implementation, plan, profile_result, fd_result,
scout_result, remaining_result (path/SHA256 pins), optional pins[].
Report additionally consumes run_result and performs zero molecular calls.
Run uses one warm native MACE model and two single-thread GFN2 workers;
at most four new composite cells. No search, force, chemistry or score changes.
"""
from concurrent.futures import ThreadPoolExecutor
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import socket
import sys
import threading
import time
import traceback

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
PROTOCOL = 'lbt3_whole_terminal_common_pool_v1'
SEARCH_PROTOCOL = 'lbt3_whole_own_composite_bounded_accommodation_v1'
SOURCES, METALS = ('7CCO', '7CCN'), ('La', 'Dy')
EV_TO_KCAL, HA_TO_KCAL = 23.06054783061903, 627.509474
DESIGN_KEYS = ('plan', 'profile_result', 'fd_result', 'scout_result', 'remaining_result')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve(); h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return {'path': str(path), 'sha256': h.hexdigest()}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed input: ' + record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False); handle.write('\n')


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def difference(a, b):
    mace = math.fsum([a['MACE_eV'], -b['MACE_eV']]) * EV_TO_KCAL
    vac = math.fsum([a['GFN2_vacuum_hartree'], -b['GFN2_vacuum_hartree']]) * HA_TO_KCAL
    water = math.fsum([a['GFN2_water_hartree'], -b['GFN2_water_hartree']]) * HA_TO_KCAL
    return {'MACE_kcal_mol': mace, 'GFN2_vacuum_kcal_mol': vac,
            'GFN2_water_kcal_mol': water, 'solvent_kcal_mol': math.fsum([water, -vac]),
            'composite_kcal_mol': math.fsum([mace, water, -vac])}


def recipe(world, source, metal):
    return {'protocol': SEARCH_PROTOCOL, 'preparation': world['preparation_pin'],
        'model': world['mcfg']['model'], 'MACE_versions': world['mcfg']['expected_versions'],
        'GFN2_runtime': world['scfg']['runtime'], 'GFN2_settings': world['manifest']['settings'],
        'GFN2_tolerances': world['scfg']['tolerances'], 'physical_charge': -2,
        'physical_multiplicity': 1 if metal == 'La' else 6, 'GFN2_effective_multiplicity': 1,
        'physical_ids': [a['id'] for a in source['atoms']],
        'symbols': [metal if a['element'] == 'M' else a['element'] for a in source['atoms']]}


def check_coordinates(world, source, metal, xyz_pin, physical_pin, positions_pin=None):
    symbols, positions = world['fd'].xyz(xyz_pin)
    atoms = read(physical_pin); expected = source['atoms']
    require(len(atoms) == 256 and symbols == recipe(world, source, metal)['symbols'], 'Candidate symbols/order differ')
    require([{k: v for k, v in a.items() if k != 'xyz_A'} for a in atoms] ==
            [{k: v for k, v in a.items() if k != 'xyz_A'} for a in expected], 'Candidate physical atom identity changed')
    require(positions == [a['xyz_A'] for a in atoms], 'XYZ/physical coordinates differ')
    x = np.asarray(positions, dtype=float)
    require(x.shape == (256, 3) and np.isfinite(x).all(), 'Candidate coordinates invalid')
    if positions_pin is not None:
        require(np.array_equal(x, np.load(verify(positions_pin), allow_pickle=False)), 'Saved physical positions differ')
    return x


def terminal_guard(world, source, search, x):
    """Replay source-relative guards, including mobile metal against fixed peptide."""
    from scipy.spatial import cKDTree
    from scaffold_accommodation import Geometry, SETTINGS
    active = world['profile']['active_candidate_basis']['mode_indices_selection_order']
    kin = source['kin']; q = np.asarray(search['full_q']); original = source['x']
    require(q.shape == (64,) and np.isfinite(q).all() and search['active_indices'] == active,
            'Terminal coordinate scope differs')
    require([kin.modes[k]['name'] for k in active] == search['active_names'] ==
            world['profile']['active_candidate_basis']['mode_names_selection_order'], 'Terminal mode names differ')
    require(np.count_nonzero(q[[k for k in range(64) if k not in active]]) == 0, 'Inactive coordinate changed')
    replay, _, _, _ = kin.evaluate(q)
    require(np.allclose(replay, x, rtol=0, atol=1e-12), 'Terminal kinematic replay differs')
    moving = set(i for k in active for i in kin.modes[k]['moving_indices']); fixed = sorted(set(range(256)) - moving)
    bonds = [b['native_particle_indices'] for b in source['bonds']]
    adapted = [dict(source_id=a['id'], residue_id=a['id'].rsplit('/', 1)[0], name=a['id'].rsplit('/', 1)[1],
                    element=a['element'], resname=a['resname']) for a in source['atoms'][:255]]
    guard = Geometry(adapted, bonds, original[:255], original[:255], sorted(moving - {255}),
                     [{'element': search['metal'], 'xyz_A': original[255].tolist()}])
    guard.other = x[255:].copy()  # guard.all_source still contains the ORIGINAL metal.
    inherited = guard.check(x[:255]); graph = [set() for _ in range(256)]; excluded = set()
    for a, b in bonds:
        graph[a].add(b); graph[b].add(a)
    for a, neighbors in enumerate(graph):
        near = set(neighbors)
        for b in neighbors:
            near.update(graph[b])
        excluded.update(tuple(sorted((a, b))) for b in near if a != b)
    overlaps = []
    for a, b in sorted(cKDTree(x).query_pairs(SETTINGS['new_overlap_heavy_A'])):
        if (a not in moving and b not in moving) or (a, b) in excluded:
            continue
        cutoff = SETTINGS['new_overlap_H_A'] if 'H' in (source['atoms'][a]['element'], source['atoms'][b]['element']) else SETTINGS['new_overlap_heavy_A']
        if np.linalg.norm(x[a] - x[b]) < min(cutoff, np.linalg.norm(original[a] - original[b])) - SETTINGS['overlap_roundoff_A']:
            overlaps.append([a, b])
    angles = np.asarray([(a, b, c) for b, ns in enumerate(graph) for a in sorted(ns) for c in sorted(ns) if a < c])
    def cosines(y):
        a, b, c = angles.T; u, v = y[a] - y[b], y[c] - y[b]
        return np.einsum('ij,ij->i', u, v) / np.linalg.norm(u, axis=1) / np.linalg.norm(v, axis=1)
    angle = float(np.max(np.abs(cosines(x) - cosines(original))))
    protected = max(float(np.max(np.abs(np.linalg.norm(x[g['indices']][:, None] - x[g['indices']][None, :], axis=2) -
                np.linalg.norm(original[g['indices']][:, None] - original[g['indices']][None, :], axis=2))))
                for g in source['modes']['protected_groups'])
    extent = float(np.linalg.norm(x - original, axis=1).max()); target = float(np.linalg.norm(x[255] - original[255]))
    fixed_error = float(np.max(np.abs(x[fixed] - original[fixed]), initial=0))
    torsion = max(abs(float(q[k])) for k in active if kin.modes[k]['unit'] == 'radian')
    require(inherited['pass'] and not overlaps and inherited['source_bond_max_error_A'] <= 1e-9
            and angle <= 1e-10 and protected <= 1e-9 and fixed_error <= 1e-12
            and extent <= .20 + 1e-12 and target <= .20 + 1e-12 and torsion <= .20 + 1e-12,
            'Terminal geometry fails original-source physical admission')
    return {'pass': True, 'baseline': source['prepared']['physical_atoms'], 'all256_contacts': True,
            'maximum_source_displacement_A': extent, 'target_displacement_A': target,
            'fixed_atom_change_A': fixed_error, 'maximum_torsion_radian': torsion,
            'bond_angle_cosine_change': angle, 'protected_distance_change_A': protected,
            'inherited_peptide_guard': inherited, 'new_worsened_contacts': overlaps}


def point_components(world, source, metal, point):
    require(point['status'] == 'complete' and point['source_id'] == source['id'] and point['metal'] == metal
            and point['recipe'] == recipe(world, source, metal), 'Own-state point recipe differs')
    x = check_coordinates(world, source, metal, point['xyz'], point['physical_source'], point['positions'])
    require(point['cache_key'] == digest({'recipe': point['recipe'], 'positions_A': x.tolist()}), 'Point cache identity differs')
    verify(point['gradient_kcal_mol_A'])
    base = source['prepared']['endpoints']['whole'][metal]
    if 'reused_origin' in point:
        require(np.array_equal(x, source['x']) and point['reused_origin'] ==
                {'MACE': world['profile']['mace_forces'], 'solvent': world['profile']['solvent_result'],
                 'origin_task_id': base['task_id']}, 'Reused origin lineage differs')
        components = world['origin_components'][(source['id'], metal)]
    else:
        mace = read(point['MACE']); state = mace['input_state_check']
        require(state['charge'] == -2 and state['spin_multiplicity'] == (1 if metal == 'La' else 6)
                and state['atoms'] == 256 and state['selected_metal_index'] == 255 and state['head'] == 'omol',
                'Own-state MACE state differs')
        components = {'MACE_eV': mace['energy_eV']}
        for medium, key in (('vacuum', 'GFN2_vacuum_hartree'), ('ALPB(Water)', 'GFN2_water_hartree')):
            low = read(point['GFN2'][medium]); task = low['task']
            require(low['status'] == 'accepted' and task['origin_task_id'] == base['task_id'] and task['medium'] == medium
                    and task['source_id'] == source['id'] and task['metal'] == metal and task['charge'] == -2
                    and task['physical_multiplicity'] == base['multiplicity'] and task['multiplicity'] == 1 and task['uhf'] == 0
                    and task['source_xyz'] == point['xyz'] and task['physical_source'] == point['physical_source'],
                    'Own-state GFN2 task differs')
            require(verify(task['xyz']).read_bytes() == verify(point['xyz']).read_bytes(), 'GFN2 scored coordinates differ')
            for record in (low['execution'], low['xtbout'], low['gradient']['artifact']):
                verify(record)
            components[key] = low['energy_hartree']
    require(components == point['components'] and all(math.isfinite(v) for v in components.values()), 'Point energy assembly differs')
    return x, components


def inputs(config_path):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    cfg = json.loads(config_path.read_text())
    require(Path(cfg['repository']).resolve() == ROOT and pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Repository/adapter differs')
    verify(cfg['plan'])
    for record in cfg.get('pins', []):
        verify(record)
    profile, fd_result, scout, remaining = [read(cfg[key]) for key in ('profile_result', 'fd_result', 'scout_result', 'remaining_result')]
    require(profile['status'] == 'complete' and fd_result['status'] == 'passed' and fd_result['composite_passed_all_states']
            and fd_result['design']['profile_result'] == cfg['profile_result'], 'Qualified profile/derivative lineage differs')
    for result in (scout, remaining):
        parent_cfg = read(result['config']); verify(result['implementation'])
        require(result['profile_result'] == cfg['profile_result'] and result['fd_result'] == cfg['fd_result']
                and result['plan'] == parent_cfg['plan'] == cfg['plan']
                and result['implementation'] == parent_cfg['implementation']
                and result['native_helper'] == fd_result['implementation'], 'Search parent/recipe lineage differs')
    require(scout['stage'] == 'scout' and remaining['stage'] == 'remaining'
            and remaining['scout_result'] == cfg['scout_result'] and remaining['settings'] == scout['settings']
            and remaining['frozen_scout_implementation'] == scout['implementation'], 'Declared search lineage differs')
    solvent = read(profile['solvent_result']); mforce = read(profile['mace_forces']); modes = read(profile['modes_result'])
    mcfg, scfg = read(mforce['config']), read(solvent['config'])
    require(solvent['status'] == mforce['status'] == modes['status'] == 'complete'
            and solvent['accepted'] == solvent['attempted'] == 8 and mforce['model'] == mcfg['model']
            and mforce['versions'] == mcfg['expected_versions']
            and solvent['preparation'] == modes['preparation'] == mcfg['preparation'], 'Origin/model lineage differs')
    manifest = read(solvent['manifest']); require(manifest['runtime'] == scfg['runtime'], 'Standalone runtime differs')
    solvent_module = module('lbt3_pool_solvent', solvent['implementation'])
    require(manifest['settings'] == solvent_module.EXPECTED_SETTINGS and scfg['tolerances'] == solvent_module.TOLERANCES,
            'Qualified standalone settings differ')
    for name, record in profile['helpers'].items():
        verify(record)
        require(pin(ROOT / 'scripts' / name)['sha256'] == record['sha256'], 'Geometry helper differs: ' + name)
    sys.path.insert(0, str(ROOT / 'scripts'))
    global np
    import numpy as np
    from mace_site_kinematics import Kinematics
    prep = read(modes['preparation']); prepared_rows = {r['source_id']: r for r in prep['rows']}
    mode_rows = {r['source_id']: r for r in modes['rows']}
    force_by = {r['task_id']: r for r in mforce['rows']}
    raw = {(r['task']['origin_task_id'], r['task']['medium']): r for r in solvent['rows']}
    require(len(force_by) == 4 and len(raw) == 8 and all(r['status'] == 'accepted' for r in raw.values()), 'Origin matrix incomplete')
    world = dict(profile=profile, fd=module('lbt3_pool_fd', fd_result['implementation']), mcfg=mcfg, scfg=scfg,
                 manifest=manifest, solvent=solvent, raw=raw, preparation_pin=modes['preparation'], origin_components={})
    sources = {}; signatures = []
    for source_id in SOURCES:
        prepared = read(prepared_rows[source_id]['preparation']); atoms = read(prepared['physical_atoms'])
        mode = read(mode_rows[source_id]['modes']); kdata = read(mode['kinematics'])
        require(prepared['status'] == 'prepared' and prepared['charge'] == -2 and prepared['protein_charge'] == -5
                and prepared['waters'] == prepared['spectators'] == prepared['caps'] == [] and len(atoms) == 256
                and mode['physical_source'] == prepared['physical_atoms'] and mode['source_preparation'] == prepared_rows[source_id]['preparation']
                and kdata['positions_A'] == [a['xyz_A'] for a in atoms] and kdata['physical_ids'] == [a['id'] for a in atoms],
                'Whole source preparation/kinematics differs')
        bonds, native = read(prepared['covalent_bonds']), read(prepared['native_peptide_atoms'])
        signatures.append({'atoms': [(a['id'], a['element'], a.get('kind'), a.get('resname')) for a in atoms],
            'charges': [(a['id'], a['charge_e']) for a in native],
            'bonds': sorted((tuple(sorted((b['atom_a_id'], b['atom_b_id']))), b['forcefield_equilibrium_A']) for b in bonds),
            'terminal_chemistry': prepared['terminal_chemistry'], 'forcefield': prepared['forcefield']})
        source = dict(id=source_id, prepared=prepared, atoms=atoms, bonds=bonds, modes=mode,
                      kin=Kinematics(kdata), x=np.asarray([a['xyz_A'] for a in atoms]))
        for metal in METALS:
            task = prepared['endpoints']['whole'][metal]; force = force_by[task['task_id']]
            require(task['source_id'] == source_id and task['charge'] == -2 and task['multiplicity'] == (1 if metal == 'La' else 6)
                    and task['atom_count'] == 256 and task['metal_index'] == 255 and task['source_cap_map'] is None
                    and task['oxidation_state'] == 3 and force['status'] == 'complete'
                    and all(force[k] == v for k, v in task.items()), 'Prepared origin/native state differs')
            check_coordinates(world, source, metal, task['xyz'], task['physical_source'])
            for medium in ('vacuum', 'ALPB(Water)'):
                low_task = raw[(task['task_id'], medium)]['task']
                require(low_task['source_id'] == source_id and low_task['metal'] == metal and low_task['charge'] == -2
                        and low_task['physical_multiplicity'] == task['multiplicity'] and low_task['multiplicity'] == 1
                        and low_task['uhf'] == 0 and low_task['source_xyz'] == task['xyz']
                        and low_task['physical_source'] == task['physical_source'], 'Origin GFN2 state/coordinates differ')
            world['origin_components'][(source_id, metal)] = {'MACE_eV': force['energy_eV'],
                'GFN2_vacuum_hartree': raw[(task['task_id'], 'vacuum')]['energy_hartree'],
                'GFN2_water_hartree': raw[(task['task_id'], 'ALPB(Water)')]['energy_hartree']}
        sources[source_id] = source
    require(signatures[0] == signatures[1], 'Across-source chemical inventory differs')
    slots = [dict(source_id=scout['source_id'], metal=scout['metal'], status=scout['status'],
                  search=scout['search'], points=scout['points'], parent=cfg['scout_result'])]
    slots.extend(dict(row, parent=cfg['remaining_result']) for row in remaining['rows'])
    require(len(slots) == 4 and {(r['source_id'], r['metal']) for r in slots} == {(s, m) for s in SOURCES for m in METALS},
            'Four-search denominator differs')
    world.update(sources=sources, slots=slots, parents={'scout': scout, 'remaining': remaining})
    return cfg, world


def assemble(world):
    candidates, by_key, aliases, cells, slots = [], {}, {}, {}, []
    def add(source_id, x, physical, alias):
        # Normalize signed zero only; no rounding, fitting or alignment.
        canonical = np.where(x == 0., 0., x).astype('<f8')
        key = hashlib.sha256(canonical.tobytes()).hexdigest()
        if key not in by_key:
            candidate = {'candidate_id': 'geometry_%02d' % len(candidates), 'coordinate_sha256': key,
                         'source_id': source_id, 'physical_source': physical, 'aliases': [], 'x': x}
            by_key[key] = candidate; candidates.append(candidate)
        candidate = by_key[key]
        require(np.array_equal(candidate['x'], x), 'Coordinate hash collision')
        candidate['aliases'].append(alias); aliases[alias] = candidate['candidate_id']
        return candidate['candidate_id']
    def reuse(ident, metal, components, evidence):
        key = (ident, metal)
        if key in cells:
            require(cells[key]['components'] == components, 'Conflicting exact geometry/state cached energies')
            cells[key]['evidence'].append(evidence)
        else:
            cells[key] = {'candidate_id': ident, 'metal': metal, 'status': 'accepted', 'components': components,
                          'reused': True, 'evidence': [evidence]}
    for source_id in SOURCES:
        source = world['sources'][source_id]
        ident = add(source_id, source['x'], source['prepared']['physical_atoms'], source_id + ':origin')
        for metal in METALS:
            reuse(ident, metal, world['origin_components'][(source_id, metal)],
                  {'origin_task': source['prepared']['endpoints']['whole'][metal],
                   'MACE': world['profile']['mace_forces'], 'solvent': world['profile']['solvent_result']})
    for slot in sorted(world['slots'], key=lambda r: (SOURCES.index(r['source_id']), METALS.index(r['metal']))):
        source_id, metal = slot['source_id'], slot['metal']; source = world['sources'][source_id]
        item = {'source_id': source_id, 'metal': metal, 'parent': slot['parent'], 'status': 'unavailable'}
        if slot.get('search') is None:
            item['failure'] = 'Search did not produce a search receipt'; slots.append(item); continue
        search = read(slot['search']); item['search'] = slot['search']
        require(search['source_id'] == source_id and search['metal'] == metal, 'Search slot identity differs')
        item.update(termination=search.get('termination'), accepted_steps=search.get('accepted_steps'))
        if slot['status'] != 'complete' or not search['candidate_admitted']:
            item['failure'] = 'Declared search unavailable; archived/rejected points are not candidates'; slots.append(item); continue
        require(search['status'] == 'finite_proposal_available' and search['final_geometry']['pass']
                and search['origin_guard']['pass'] and search['candidate'] in slot['points']
                and search['candidate'] == search['trace'][-1]['point']
                and len(search['trace']) == search['accepted_steps'] + 1
                and all(r['geometry']['pass'] for r in search['trace'])
                and search['authoritative_origin'] == source['prepared']['endpoints']['whole'][metal]['xyz'],
                'Terminal proposal admission/provenance differs')
        point = read(search['candidate']); x, components = point_components(world, source, metal, point)
        geometry = terminal_guard(world, source, search, x)
        ident = add(source_id, x, point['physical_source'], source_id + ':terminal_' + metal)
        reuse(ident, metal, components, {'terminal_point': search['candidate'], 'search': slot['search'], 'parent': slot['parent']})
        work = difference(components, world['origin_components'][(source_id, metal)])
        require(work == search['work_from_origin'] and work['composite_kcal_mol'] <= 1e-8, 'Own terminal work differs')
        item.update(status='admitted', candidate_id=ident, point=search['candidate'], physical_check=geometry,
                    own_metal_work=work, final_force_diagnostic=search['final_force_diagnostic'],
                    boundary_flags=search['preserved_boundary_flags'])
        slots.append(item)
    missing = [(c['candidate_id'], m) for c in candidates for m in METALS if (c['candidate_id'], m) not in cells]
    require(len(candidates) <= 6 and len(missing) <= 4, 'Declared terminal pool/cross-cell ceiling exceeded')
    return candidates, aliases, cells, slots, missing


def execute(world, candidates, missing, out, receipt):
    require(len(missing) <= 4 and int(os.environ['SLURM_CPUS_PER_TASK']) >= 2, 'Cross-call/allocation ceiling differs')
    fd, scfg, mcfg = world['fd'], world['scfg'], world['mcfg']
    runtime = scfg['runtime']; runtime_pins = [runtime['executable'], runtime['parameters'], *runtime['runtime_libraries']]
    for record in runtime_pins:
        verify(record)
    backend = module('lbt3_pool_backend', scfg['backend_helper']); env = os.environ.copy()
    for key, value in world['manifest']['effective_environment'].items():
        if value is None:
            env.pop(key, None)
        else:
            env[key] = value
    require(env.get('OMP_NUM_THREADS') == '1', 'Standalone thread count differs')
    native = None; rows = []; by = {c['candidate_id']: c for c in candidates}
    cost = receipt['actual_cost']; cost_lock = threading.Lock()
    for ident, metal in missing:
        candidate = by[ident]; source = world['sources'][candidate['source_id']]
        base = source['prepared']['endpoints']['whole'][metal]; directory = out / (ident + '_' + metal); directory.mkdir()
        symbols = recipe(world, source, metal)['symbols']; x = candidate['x']
        physical = read(candidate['physical_source']); require(x.tolist() == [a['xyz_A'] for a in physical], 'Cross-cell physical mapping differs')
        with (directory / 'whole.xyz').open('x') as handle:
            handle.write('256\nLBT3 common-pool fixed physical geometry\n')
            for symbol, position in zip(symbols, x):
                handle.write(symbol + ' ' + ' '.join(format(v, '.17g') for v in position) + '\n')
        xyz_pin = pin(directory / 'whole.xyz')
        task = dict(base, task_id=ident + '_' + metal, xyz=xyz_pin, physical_source=candidate['physical_source'],
                    origin_task_id=base['task_id'], candidate_id=ident)
        check_coordinates(world, source, metal, xyz_pin, candidate['physical_source'])
        write(directory / 'TASK.json', task)
        row = {'candidate_id': ident, 'metal': metal, 'status': 'unavailable', 'reused': False,
               'task': pin(directory / 'TASK.json'), 'recipe': recipe(world, source, metal),
               'cache_key': digest({'recipe': recipe(world, source, metal), 'positions_A': x.tolist()}),
               'new_MACE_calls': 0, 'new_GFN2_calls': 0}
        cost['new_composite_evaluations'] += 1; tick = time.monotonic()
        try:
            if native is None:
                load = time.monotonic(); native = fd.NativeMACE(mcfg); cost['model_load_seconds'] += time.monotonic() - load
            before = native.calls; start = time.monotonic()
            try:
                mace = native.evaluate(task, symbols, x, compute_force=False)
            finally:
                row['new_MACE_calls'] = native.calls - before; cost['MACE_calls_attempted'] += row['new_MACE_calls']
                cost['MACE_call_seconds'] += time.monotonic() - start
            write(directory / 'MACE.json', mace); row['MACE'] = pin(directory / 'MACE.json')
            def low(medium):
                target = directory / ('water' if medium == 'ALPB(Water)' else 'vacuum'); target.mkdir()
                (target / 'core.xyz').write_bytes(verify(xyz_pin).read_bytes())
                old = world['raw'][(base['task_id'], medium)]['task']
                command = [runtime['executable']['path'], 'core.xyz', '--gfn', '2', '--chrg', '-2', '--uhf', '0',
                           '--etemp', '300', '--acc', '0.0001', '--iterations', '500', '--norestart', '--json', '--pop', '--sp']
                if medium == 'ALPB(Water)':
                    command += ['--alpb', 'water', 'gsolv']
                local = dict(old, task_id=task['task_id'] + '_' + target.name, directory=str(target), command=command,
                             xyz=pin(target / 'core.xyz'), source_xyz=xyz_pin, physical_source=candidate['physical_source'], gradient_requested=False)
                local.pop('task_record', None); write(target / 'TASK.json', local); local['task_record'] = pin(target / 'TASK.json')
                result = {'status': 'unavailable', 'task': local, 'call_attempted': False}; started = time.monotonic()
                try:
                    for record in [local['task_record'], local['xyz'], *runtime_pins]:
                        verify(record)
                    result['call_attempted'] = True
                    with cost_lock:
                        cost['GFN2_calls_attempted'] += 1
                    execution = backend.call(command, target, env)
                    result['execution'] = pin(target / 'EXECUTION.json'); verify(local['xyz'])
                    result.update(fd.admit_energy(target, local, execution, backend, scfg['tolerances']), status='accepted')
                    for record in runtime_pins:
                        verify(record)
                except Exception as exc:
                    result.update(failure=repr(exc), traceback=traceback.format_exc())
                result['elapsed_seconds'] = time.monotonic() - started
                with cost_lock:
                    cost['GFN2_call_seconds'] += result['elapsed_seconds']
                result['artifacts'] = {p.name: pin(p) for p in target.iterdir() if p.is_file()}
                write(target / 'RESULT.json', result)
                return result
            before_low = cost['GFN2_calls_attempted']
            try:
                with ThreadPoolExecutor(max_workers=2) as pool:
                    low_rows = list(pool.map(low, ('vacuum', 'ALPB(Water)')))
            finally:
                row['new_GFN2_calls'] = cost['GFN2_calls_attempted'] - before_low
            row['GFN2'] = {r['task']['medium']: pin(Path(r['task']['directory']) / 'RESULT.json') for r in low_rows}
            require(all(r['status'] == 'accepted' for r in low_rows), 'Cross-state GFN2 component unavailable')
            row.update(status='accepted', components={'MACE_eV': mace['energy_eV'],
                'GFN2_vacuum_hartree': low_rows[0]['energy_hartree'], 'GFN2_water_hartree': low_rows[1]['energy_hartree']})
        except Exception as exc:
            row.update(failure=repr(exc), traceback=traceback.format_exc())
        row['elapsed_seconds'] = time.monotonic() - tick; write(directory / 'CELL.json', row)
        row['cell'] = pin(directory / 'CELL.json'); rows.append(row)
        receipt['cross_rows'] = rows
        print(json.dumps({'candidate_id': ident, 'metal': metal, 'status': row['status']}), flush=True)
    if native is not None:
        receipt['native_runtime'] = {'versions': native.versions, 'device': native.device, 'actual_model_forwards': native.calls,
            'peak_cuda_allocated_bytes': native.torch.cuda.max_memory_allocated()}
    require(cost['new_composite_evaluations'] <= 4 and cost['MACE_calls_attempted'] <= 4 and cost['GFN2_calls_attempted'] <= 8,
            'Cross-cell call ceiling exceeded')
    return rows


def pool_summary(ids, cells):
    ids = list(dict.fromkeys(ids))
    missing = [{'candidate_id': ident, 'metal': m} for ident in ids for m in METALS if (ident, m) not in cells]
    if missing or not ids:
        return {'status': 'unavailable', 'candidate_ids': ids, 'missing_cells': missing}
    selected, ties, work = {}, {}, {}
    for metal in METALS:
        reference = cells[(ids[0], metal)]['components']
        values = {ident: difference(cells[(ident, metal)]['components'], reference)['composite_kcal_mol'] for ident in ids}
        best = min(values.values()); ties[metal] = [ident for ident in ids if values[ident] == best]
        selected[metal] = ties[metal][0]; work[metal] = values
    paired = {ident: difference(cells[(ident, 'Dy')]['components'], cells[(ident, 'La')]['components']) for ident in ids}
    return {'status': 'available', 'candidate_ids': ids, 'selected': selected, 'exact_ties': ties,
        'selection_rule': 'each metal minimizes same complete composite menu; exact ties use first declared candidate',
        'same_metal_work_relative_first_kcal_mol': work, 'fixed_geometry_Dy_minus_La': paired,
        'independently_selected_Dy_minus_La': difference(cells[(selected['Dy'], 'Dy')]['components'], cells[(selected['La'], 'La')]['components']),
        'La_selected_geometry_Dy_minus_La': paired[selected['La']], 'Dy_selected_geometry_Dy_minus_La': paired[selected['Dy']]}


def summarize(candidates, aliases, cells, slots):
    original_ids = [aliases[s + ':origin'] for s in SOURCES]
    result = {'original_two_source_common_pool': pool_summary(original_ids, cells), 'source_local_pools': {},
              'contrast_convention': 'Dy minus La; comparison of these contrasts across geometries is conditional energy response, not affinity',
              'source_spread_convention': '7CCN minus 7CCO and absolute difference of source-local independently selected contrasts',
              'global_pool_source_spread_not_claimed': 'One shared global menu removes source choice by construction; source-local menus test remaining source sensitivity.'}
    for source in SOURCES:
        origin = aliases[source + ':origin']; local = {'origin': pool_summary([origin], cells)}
        for label, metals in (('plus_La_terminal', ('La',)), ('plus_Dy_terminal', ('Dy',)), ('plus_both_terminals', METALS)):
            names = [source + ':terminal_' + metal for metal in metals]
            local[label] = pool_summary([origin] + [aliases[n] for n in names], cells) if all(n in aliases for n in names) else {
                'status': 'unavailable', 'reason': 'A declared terminal search is unavailable'}
        result['source_local_pools'][source] = local
    result['source_spread'] = {}
    for menu in ('origin', 'plus_La_terminal', 'plus_Dy_terminal', 'plus_both_terminals'):
        pair = [result['source_local_pools'][s][menu] for s in SOURCES]
        if all(p['status'] == 'available' for p in pair):
            a, b = [p['independently_selected_Dy_minus_La']['composite_kcal_mol'] for p in pair]
            result['source_spread'][menu] = {'7CCN_minus_7CCO_kcal_mol': math.fsum([b, -a]), 'absolute_spread_kcal_mol': abs(b-a)}
        else:
            result['source_spread'][menu] = {'status': 'unavailable'}
    all_searches = all(s['status'] == 'admitted' for s in slots)
    result['expanded_global_common_pool'] = pool_summary([c['candidate_id'] for c in candidates], cells) if all_searches else {
        'status': 'unavailable', 'reason': 'Four declared terminal searches are required; original pool preserved separately'}
    original, expanded = result['original_two_source_common_pool'], result['expanded_global_common_pool']
    if expanded['status'] == 'available':
        result['global_pool_change'] = {'selection_changed': {m: original['selected'][m] != expanded['selected'][m] for m in METALS},
            'own_metal_work': {m: difference(cells[(expanded['selected'][m], m)]['components'], cells[(original['selected'][m], m)]['components']) for m in METALS},
            'Dy_minus_La_shift_kcal_mol': math.fsum([expanded['independently_selected_Dy_minus_La']['composite_kcal_mol'],
                                                   -original['independently_selected_Dy_minus_La']['composite_kcal_mol']])}
    if 'absolute_spread_kcal_mol' in result['source_spread']['plus_both_terminals']:
        result['source_spread_change_kcal_mol'] = result['source_spread']['plus_both_terminals']['absolute_spread_kcal_mol'] - result['source_spread']['origin']['absolute_spread_kcal_mol']
    result['paired_matrix'] = [{'candidate_id': c['candidate_id'], 'aliases': c['aliases'],
        'cells': {m: cells.get((c['candidate_id'], m), {'status': 'unavailable'}) for m in METALS}} for c in candidates]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('run', 'report'))
    for name in ('config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args(); cfg, world = inputs(args.config.resolve())
    require(pin(sys.executable)['sha256'] == world['mcfg' if args.stage == 'run' else 'scfg']['python']['sha256'], 'Stage Python differs')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    design = {k: cfg[k] for k in DESIGN_KEYS}; design['implementation_sha256'] = cfg['implementation']['sha256']
    receipt = {'protocol_id': PROTOCOL, 'stage': args.stage, 'status': 'unavailable', 'design': design,
        'implementation': pin(__file__), 'config': pin(args.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
        'host': socket.gethostname(), 'cross_rows': [], 'new_gradients_requested': 0, 'optimizer_steps': 0,
        'actual_cost': {'new_composite_evaluations': 0, 'MACE_calls_attempted': 0, 'GFN2_calls_attempted': 0,
                        'model_load_seconds': 0., 'MACE_call_seconds': 0., 'GFN2_call_seconds': 0.},
        'interpretation_limits': 'Fixed inventory and conditional termini; bounded local response. No binding free energy, Kd, equilibrium, full stationarity or Mex composite claim.'}
    write(out / 'START.json', receipt)
    try:
        candidates, aliases, cells, slots, missing = assemble(world)
        receipt.update(candidates=[{k: v for k, v in c.items() if k != 'x'} for c in candidates], aliases=aliases,
            searches=slots, declared_searches=4, declared_slots=6, required_cross_cells=[{'candidate_id': i, 'metal': m} for i, m in missing],
            original_summary=pool_summary([aliases[s + ':origin'] for s in SOURCES], cells))
        write(out / 'POOL.json', {k: receipt[k] for k in ('design', 'candidates', 'aliases', 'searches', 'required_cross_cells')})
        receipt['pool'] = pin(out / 'POOL.json')
        if args.stage == 'run':
            rows = execute(world, candidates, missing, out, receipt) if all(s['status'] == 'admitted' for s in slots) else []
        else:
            previous = read(cfg['run_result']); read(previous['config']); verify(previous['implementation'])
            require(previous['stage'] == 'run' and previous['status'] in ('complete', 'unavailable') and previous['design'] == design
                    and previous['candidates'] == receipt['candidates'] and previous['aliases'] == aliases
                    and previous['required_cross_cells'] == receipt['required_cross_cells'], 'Run/pool design differs')
            rows = previous['cross_rows']; receipt['run_result'] = cfg['run_result']
        require(len(rows) == len({(r['candidate_id'], r['metal']) for r in rows}) and
                {(r['candidate_id'], r['metal']) for r in rows} <= set(missing), 'Unexpected/duplicate cross cells')
        for row in rows:
            require(read(row['cell']) == {k: v for k, v in row.items() if k != 'cell'}, 'Cross-cell/parent row differs')
            if row['status'] == 'accepted':
                require(all(math.isfinite(v) for v in row['components'].values()), 'Nonfinite cross energy')
                cells[(row['candidate_id'], row['metal'])] = row
        receipt['cross_rows'] = rows; receipt['summary'] = summarize(candidates, aliases, cells, slots)
        receipt['status'] = 'complete' if receipt['summary']['expanded_global_common_pool']['status'] == 'available' else 'unavailable'
        receipt['parent_call_costs'] = {k: p['actual_cost'] for k, p in world['parents'].items()}
        verify(receipt['config'])
    except Exception as exc:
        receipt.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    receipt['elapsed_seconds'] = time.monotonic() - started
    receipt['application_allocated_CPU_seconds'] = receipt['elapsed_seconds'] * int(os.environ['SLURM_CPUS_PER_TASK'])
    receipt['application_allocated_GPU_seconds'] = receipt['elapsed_seconds'] if args.stage == 'run' else 0
    receipt['peak_host_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    receipt['allocation'] = {k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_CPU', 'CUDA_VISIBLE_DEVICES')}
    write(out / 'RESULT.json', receipt)
    print(json.dumps({'status': receipt['status'], 'result': pin(out / 'RESULT.json')}), flush=True)
    return 0 if receipt['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
