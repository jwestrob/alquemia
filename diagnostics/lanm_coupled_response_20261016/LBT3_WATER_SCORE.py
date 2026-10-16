#!/usr/bin/env python3
"""Finite 259-atom LBT3 water energies; pinned native engines, Slurm only."""
from __future__ import annotations
import argparse
import ast
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import inspect
import json
import math
import os
from pathlib import Path
import resource
import sys
import textwrap
import threading
import time
import traceback

PROTOCOL = 'lbt3_fixed_source_one_water_composite_v1'
PLAN_SHA = '1eda4edc5bd5da4103b88e11ee75aff0a35dfeb00bf91f1362996743cb14bd10'
PREPARE_SHA = '4da7ab6c9a1600919155ed1079bae3a282de13150f8d806596a5b1d6ecbf30d1'
FD_SHA = '0e48b8c6c2f271512f2c445f826036d8a68f7dca0aabc5923697ab3ecb7d9896'
SOLVENT_SHA = '517097f1db17fc7ed1067265f3af91b6203c52c2b91eeaf62e4cc47c38e0fd90'
FD_RESULT_SHA = 'c1b88ea94c0eb9ecfcad939f13fa89d2ae8356fb859b6cf4e4cf3a37bc6117fa'
SOURCES, METALS = ('7CCO', '7CCN'), ('La', 'Dy')
POSES = ('inner_plane0', 'inner_plane90', 'outer_plane0', 'outer_plane90')
SCOUT = ('7CCO_inner_plane0_La', '7CCO_inner_plane0_Dy')
EV, HA = 23.06054783061903, 627.509474
DESIGN = ('plan', 'water_preparation_result', 'solvent_result', 'fd_result')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve(); digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return {'path': str(path), 'sha256': digest.hexdigest()}


def verify(value):
    path = Path(value['path']).resolve(strict=True)
    require(path.is_file() and pin(path)['sha256'] == value['sha256'], 'changed input: ' + value['path'])
    return path


def read(value):
    return json.loads(verify(value).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False); stream.write('\n')


def module(name, value):
    spec = importlib.util.spec_from_file_location(name, verify(value))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def atom_count_adapter(function, expected_count):
    """Copy a pinned function; change only declared integer atom-count nodes."""
    source = textwrap.dedent(inspect.getsource(function)); tree = ast.parse(source)
    nodes = [node for node in ast.walk(tree) if isinstance(node, ast.Constant)
             and type(node.value) is int and node.value == 256]
    require(len(nodes) == expected_count, 'pinned function atom-count assumptions changed: ' + function.__qualname__)
    before = ast.dump(tree, include_attributes=False)
    for node in nodes:
        node.value = 259
    after = ast.dump(tree, include_attributes=False)
    namespace = dict(function.__globals__)
    exec(compile(ast.fix_missing_locations(tree), '<private_259_atom_adapter>', 'exec'), namespace)
    return namespace[function.__name__], {'function': function.__qualname__, 'old_atom_count': 256,
        'new_atom_count': 259, 'integer_nodes_replaced': expected_count,
        'original_AST_sha256': hashlib.sha256(before.encode()).hexdigest(),
        'adapted_AST_sha256': hashlib.sha256(after.encode()).hexdigest(),
        'private_copy_only': True, 'frozen_source_or_module_mutated': False}


def coordinates(task, prep_module):
    rows = prep_module.read_xyz(verify(task['xyz']))
    require(len(rows) == 259, '259-atom coordinates required')
    symbols = [r[0] for r in rows]; x = [list(r[1:]) for r in rows]
    physical = read(task['physical_source'])
    require(len(physical) == 259 and symbols ==
            [task['metal'] if a['element'] == 'M' else a['element'] for a in physical]
            and x == [a['xyz_A'] for a in physical], 'task XYZ/physical map differs')
    require(symbols[255] == task['metal'] and symbols[256:] == ['O', 'H', 'H'] and
            symbols.count(task['metal']) == 1 and all(math.isfinite(v) for row in x for v in row),
            'physical species/order/nonfinite coordinates')
    return symbols, x


def inputs(config, stage):
    import numpy as np
    cfg = json.loads(Path(config).read_text())
    require(cfg['implementation']['sha256'] == pin(__file__)['sha256'], 'score implementation differs')
    verify(cfg['implementation'])
    for key, expected in (('plan', PLAN_SHA), ('solvent_result', SOLVENT_SHA), ('fd_result', FD_RESULT_SHA)):
        require(cfg[key]['sha256'] == expected, 'declared lineage differs: ' + key); verify(cfg[key])
    for value in cfg.get('pins', []):
        verify(value)
    repo = Path(cfg['repository']).resolve(strict=True)
    sys.path.insert(0, str(repo / 'scripts'))
    water = read(cfg['water_preparation_result']); water_cfg = read(water['config'])
    require(water['status'] == 'complete' and water['protocol_id'] == 'lbt3_fixed_source_one_water_geometry_v1'
            and water['implementation']['sha256'] == PREPARE_SHA and
            water['implementation'] == water_cfg['implementation'] and water['plan'] == water_cfg['plan'],
            'declared water preparation unavailable')
    wp = module('lbt3_water_preparation_pinned', water['implementation'])
    _, originals, template, WaterCoordinates = wp.inputs(verify(water['config']))
    require(water['settings'] == wp.SETTINGS and [r['source_id'] for r in water['rows']] == list(SOURCES),
            'water geometry declaration differs')
    require(Path(water_cfg['repository']).resolve() == repo, 'water repository differs')
    by_original = {r['preparation']['source_id']: r for r in originals}
    tasks = water['tasks']; expected = [s + '_' + p + '_' + m for s in SOURCES for p in POSES for m in METALS]
    require([t['task_id'] for t in tasks] == expected and len(set(expected)) == 16, 'full sixteen-slot denominator differs')
    data = {}; preflight = []
    for row in water['rows']:
        saved = read(row['result'])
        require(saved == {k: v for k, v in row.items() if k != 'result'}, 'water source row/receipt differs')
        source = by_original[row['source_id']]; original = source['atoms']; parent = source['preparation']
        require(row['source_preparation'] == source['receipt']['preparation'] and
                row['original_physical_source'] == parent['physical_atoms'], 'water original source binding differs')
        placement = read(row['placement']); candidates = placement['candidates']
        require(len(candidates) == 24 and row['selected_candidate_index'] == placement['selected_candidate_index'],
                'water placement catalogue differs')
        selected = max(candidates, key=lambda r: r['minimum_normalized_clearance'])
        require(selected['candidate_index'] == row['selected_candidate_index'], 'geometry-only selection differs')
        mapping = read(row['source_map']); graph = read(row['graph'])
        require(len(mapping) == 259 and [r['source_atom'] for r in mapping[:256]] == original
                and graph['whole_peptide_bonds'] == source['bonds']
                and graph['water_bonds'] == [[256, 257], [256, 258]] and graph['metal_covalent_bonds'] == [],
                'fixed original graph/water map differs')
        helper = WaterCoordinates([(a['element'], *a['xyz_A']) for a in original] + template,
                                  [{'role': 'variable', 'indices': [256, 257, 258]}])
        for pose in row['poses']:
            require(read(pose['result']) == {k: v for k, v in pose.items() if k != 'result'}, 'pose receipt differs')
            physical = read(pose['physical_source']); y = np.asarray([a['xyz_A'] for a in physical])
            require(physical[:256] == original and y.shape == (259, 3) and np.isfinite(y).all(),
                    'original256 source coordinates/atoms changed')
            replay = helper.positions(pose['rigid_coordinate_q'])
            require(np.allclose(replay, y, atol=1e-12, rtol=0), 'rigid water geometry replay differs')
            contact = wp.contacts(original, y)
            require(contact['pass'] == pose['admitted'] == pose['contact_check']['pass'] and
                    abs(contact['minimum_normalized_clearance'] - pose['contact_check']['minimum_normalized_clearance']) < 1e-10,
                    'water geometry contact admission differs')
            pair = [next(t for t in tasks if t['source_id'] == row['source_id'] and
                         t['pose_id'] == pose['pose_id'] and t['metal'] == m) for m in METALS]
            for task in pair:
                require(task in row['tasks'] and task['admitted'] == pose['admitted'] and task['charge'] == -2
                        and task['multiplicity'] == task['physical_multiplicity'] == task['spin_multiplicity'] ==
                        (1 if task['metal'] == 'La' else 6) and task['metal_index'] == 255 and task['atom_count'] == 259
                        and task['water_indices'] == [256, 257, 258] and task['water_count'] == 1
                        and task['physical_source'] == pose['physical_source'] and task['source_map'] == row['source_map']
                        and task['graph'] == row['graph'] and task['placement'] == row['placement']
                        and task['original_endpoint'] == parent['endpoints']['whole'][task['metal']],
                        'paired water state/source binding differs')
                data[task['task_id']] = coordinates(task, wp)
            a, b = [data[t['task_id']] for t in pair]
            require(a[1] == b[1] and a[0][:255] == b[0][:255] and a[0][256:] == b[0][256:],
                    'paired water coordinates differ')
            preflight.append({'source_id': row['source_id'], 'pose_id': pose['pose_id'],
                              'admitted': pose['admitted'], 'exact_original256': True, 'water_replay_passed': True})
    solvent = read(cfg['solvent_result']); fd_result = read(cfg['fd_result'])
    require(solvent['status'] == 'complete' and solvent['accepted'] == solvent['attempted'] == 8 and
            fd_result['status'] == 'passed' and fd_result['composite_passed_all_states'] and
            fd_result['design']['solvent_result'] == cfg['solvent_result'], 'qualified energy lineage unavailable')
    scfg = read(solvent['config']); manifest = read(solvent['manifest'])
    force = read(solvent['mace_forces']); mcfg = read(force['config'])
    require(force['status'] == 'complete' and force['model'] == mcfg['model'] and
            force['versions'] == mcfg['expected_versions'] and scfg['preparation'] == water['preparation_result']
            and mcfg['preparation'] == water['preparation_result'], 'original model/preparation lineage differs')
    for obj in (solvent, fd_result):
        verify(obj['implementation']); verify(obj['config'])
    # Actual12500 records its implementation in CONFIG, not its RESULT.
    verify(mcfg['implementation'])
    sm = module('lbt3_water_solvent_pinned', solvent['implementation'])
    require(manifest['runtime'] == scfg['runtime'] and manifest['settings'] == sm.EXPECTED_SETTINGS
            and scfg['tolerances'] == sm.TOLERANCES, 'qualified runtime/settings/tolerances differ')
    require(fd_result['implementation']['sha256'] == FD_SHA, 'native/helper implementation differs')
    fd = module('lbt3_water_fd_pinned', fd_result['implementation'])
    evaluator, native_adaptation = atom_count_adapter(fd.NativeMACE.evaluate, 2)
    admission, low_adaptation = atom_count_adapter(fd.admit_energy, 1)
    Native = type('WaterNativeMACE259', (fd.NativeMACE,), {'evaluate': evaluator})
    qualification = read(scfg['qualification']); qm = read(qualification['manifest']); profile = read(qm['profile'])
    require(qualification['status'] == 'passed' and qualification['accepted'] == qualification['attempted'] == 8,
            'actual parameter qualification unavailable')
    verify(scfg['backend_helper']); backend = module('lbt3_water_backend_pinned', scfg['backend_helper'])
    common = module('lbt3_water_common_pinned', scfg['common_helper'])
    require(common.HA_TO_KCAL == HA and sm.EV_TO_KCAL == fd.EV_TO_KCAL == EV, 'energy units differ')
    runtime = scfg['runtime']; runtime_pins = [runtime['executable'], runtime['parameters'], *runtime['runtime_libraries']]
    for value in runtime_pins:
        verify(value)
    for value in mcfg['helper_pins'].values():
        verify(value)
    verify(mcfg['model']['checkpoint']); verify(mcfg['python'])
    require(pin(sys.executable)['sha256'] == mcfg['python']['sha256'], 'native Python differs')
    env = os.environ.copy()
    for key, value in manifest['effective_environment'].items():
        if value is None:
            env.pop(key, None)
        else:
            env[key] = value
    require(all(env.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')),
            'standalone single-thread contract differs')
    for task in tasks:
        symbols = data[task['task_id']][0]
        electrons = math.fsum(math.fsum(profile['element'][s]['refocc']) for s in symbols) + 2
        require(electrons == int(electrons) and int(electrons) % 2 == 0, 'GFN2 effective electron parity differs')
        data[task['task_id']] = (*data[task['task_id']], int(electrons))
    world = dict(cfg=cfg, tasks=tasks, coordinates=data, scfg=scfg, mcfg=mcfg, manifest=manifest,
                 Native=Native, admit=admission, backend=backend, env=env, runtime_pins=runtime_pins,
                 preflight=preflight, preparation_module=wp,
                 adaptations=[native_adaptation, low_adaptation], native_helper=fd_result['implementation'])
    return world


def component_difference(a, b):
    mace = math.fsum([a['MACE_eV'], -b['MACE_eV']]) * EV
    vacuum = math.fsum([a['GFN2_vacuum_hartree'], -b['GFN2_vacuum_hartree']]) * HA
    water = math.fsum([a['GFN2_water_hartree'], -b['GFN2_water_hartree']]) * HA
    return {'MACE_kcal_mol': mace, 'GFN2_vacuum_kcal_mol': vacuum,
            'GFN2_water_kcal_mol': water, 'solvent_kcal_mol': math.fsum([water, -vacuum]),
            'composite_kcal_mol': math.fsum([mace, water, -vacuum])}


def validate_saved_cell(world, row, task):
    require(row['status'] == 'complete' and row['task'] == task and task['admitted'], 'scout cell unavailable/state differs')
    require(read(row['cell']) == {k: v for k, v in row.items() if k != 'cell'}, 'scout CELL/row differs')
    mace = read(row['MACE']); state = mace['input_state_check']
    require(state['atoms'] == 259 and state['charge'] == -2 and state['spin_multiplicity'] == task['multiplicity']
            and state['selected_metal_index'] == 255 and state['head'] == 'omol'
            and mace['forces_eV_A'] is None and math.isfinite(mace['energy_eV']), 'scout native state differs')
    values = {'MACE_eV': mace['energy_eV']}
    for medium, key in (('vacuum', 'GFN2_vacuum_hartree'), ('ALPB(Water)', 'GFN2_water_hartree')):
        low = read(row['GFN2'][medium]); local = low['task']; directory = Path(local['directory'])
        require(low['status'] == 'accepted' and local['water_task_id'] == task['task_id']
                and local['source_xyz'] == task['xyz'] and local['physical_source'] == task['physical_source']
                and local['physical_multiplicity'] == task['multiplicity'] and local['multiplicity'] == 1
                and local['spin_multiplicity'] == 1 and local['uhf'] == 0 and local['charge'] == -2
                and local['atom_count'] == 259 and local['metal'] == task['metal'] and local['medium'] == medium
                and local['valence_electrons'] == world['coordinates'][task['task_id']][2], 'scout effective state differs')
        require(read(local['task_record']) == {k: v for k, v in local.items() if k != 'task_record'}, 'GFN2 TASK pin differs')
        require(verify(local['xyz']).read_bytes() == verify(task['xyz']).read_bytes(), 'scout scored coordinates differ')
        for value in low['artifacts'].values():
            verify(value)
        execution = read(low['execution'])
        actual = world['admit'](directory, local, execution, world['backend'], world['scfg']['tolerances'])
        require(actual['energy_hartree'] == low['energy_hartree'] and actual['xtbout'] == low['xtbout'],
                'scout raw energy admission replay differs')
        values[key] = actual['energy_hartree']
    require(values == row['components'], 'scout composite component binding differs')


def previous_scout(world):
    cfg = world['cfg']; prior = read(cfg['scout_result']); pcfg = read(prior['config'])
    require(prior['stage'] == 'scout' and prior['status'] == 'complete' and prior['scout_gate_passed']
            and prior['implementation'] == pcfg['implementation'] and
            prior['implementation']['sha256'] == cfg['implementation']['sha256']
            and all(prior[k] == pcfg[k] == cfg[k] for k in DESIGN), 'exact admitted paired scout required')
    verify(prior['implementation'])
    require(prior['model'] == world['mcfg']['model'] and prior['native_runtime']['versions'] == world['mcfg']['expected_versions']
            and prior['GFN2_runtime'] == world['scfg']['runtime'] and prior['GFN2_settings'] == world['manifest']['settings']
            and prior['private_atom_count_adaptations'] == world['adaptations'], 'scout physical recipe differs')
    rows = prior['rows']; require([r['task']['task_id'] for r in rows] == [t['task_id'] for t in world['tasks']],
                                 'scout sixteen-slot denominator differs')
    complete = {r['task']['task_id']: r for r in rows if r['status'] == 'complete'}
    require(set(complete) == set(SCOUT), 'scout executed scope differs')
    for task in world['tasks']:
        if task['task_id'] in complete:
            validate_saved_cell(world, complete[task['task_id']], task)
    require(prior['actual_cost']['MACE_calls_attempted'] == 2 and prior['actual_cost']['GFN2_calls_attempted'] == 4,
            'scout actual call count differs')
    return prior, complete


def execute(world, tasks, out, receipt):
    runtime = world['scfg']['runtime']; native = None; stop = None; rows = []
    costs = receipt['actual_cost']; lock = threading.Lock()
    for task in tasks:
        row = {'task': task, 'status': 'unavailable', 'new_MACE_calls': 0, 'new_GFN2_calls': 0}
        if not task['admitted']:
            row.update(status='geometry_rejected', reason='Frozen preparation rejected this slot; no replacement')
            rows.append(row); continue
        if stop:
            row['reason'] = 'Not attempted after native/process failure: ' + stop
            rows.append(row); continue
        directory = out / task['task_id']; directory.mkdir(); tick = time.monotonic()
        try:
            symbols, x, electrons = world['coordinates'][task['task_id']]
            verify(task['xyz']); write(directory / 'TASK.json', task)
            try:
                if native is None:
                    began = time.monotonic()
                    try:
                        native = world['Native'](world['mcfg'])
                    finally:
                        costs['model_load_seconds'] += time.monotonic() - began
                before = native.calls; began = time.monotonic()
                try:
                    mace = native.evaluate(task, symbols, x, compute_force=False)
                finally:
                    row['new_MACE_calls'] = native.calls - before
                    costs['MACE_calls_attempted'] += row['new_MACE_calls']
                    costs['MACE_call_seconds'] += time.monotonic() - began
                require(mace['forces_eV_A'] is None, 'energy-only native call returned forces')
            except Exception as exc:
                stop = repr(exc); raise
            write(directory / 'MACE.json', mace); row['MACE'] = pin(directory / 'MACE.json')
            verify(task['xyz'])

            def low(medium):
                target = directory / ('water' if medium == 'ALPB(Water)' else 'vacuum'); target.mkdir()
                (target / 'core.xyz').write_bytes(verify(task['xyz']).read_bytes())
                command = [runtime['executable']['path'], 'core.xyz', '--gfn', '2', '--chrg', '-2', '--uhf', '0',
                           '--etemp', '300', '--acc', '0.0001', '--iterations', '500', '--norestart', '--json', '--pop', '--sp']
                if medium == 'ALPB(Water)':
                    command += ['--alpb', 'water', 'gsolv']
                local = dict(task, task_id=task['task_id'] + '_' + target.name, water_task_id=task['task_id'],
                             directory=str(target), command=command, xyz=pin(target / 'core.xyz'), source_xyz=task['xyz'],
                             physical_multiplicity=task['multiplicity'], multiplicity=1, spin_multiplicity=1,
                             uhf=0, valence_electrons=electrons, medium=medium, gradient_requested=False,
                             physical_electronic_state=task['original_endpoint']['electronic_state'],
                             electronic_state='GFN2 spin-independent f-in-core effective singlet')
                write(target / 'TASK.json', local); local['task_record'] = pin(target / 'TASK.json')
                result = {'status': 'unavailable', 'task': local, 'call_attempted': False}
                began = time.monotonic()
                try:
                    for value in [local['task_record'], local['xyz'], *world['runtime_pins']]:
                        verify(value)
                    result['call_attempted'] = True
                    with lock:
                        costs['GFN2_calls_attempted'] += 1
                    execution = world['backend'].call(command, target, world['env'])
                    result['execution'] = pin(target / 'EXECUTION.json')
                    result.update(world['admit'](target, local, execution, world['backend'], world['scfg']['tolerances']),
                                  status='accepted')
                    verify(local['xyz'])
                    for value in world['runtime_pins']:
                        verify(value)
                except Exception as exc:
                    result.update(reason=repr(exc), traceback=traceback.format_exc())
                if (target / 'EXECUTION.json').exists():
                    result['execution'] = pin(target / 'EXECUTION.json'); execution = read(result['execution'])
                    result['fatal_process_state'] = not execution['owned_process_group_empty']
                    with lock:
                        costs['GFN2_solver_wall_seconds'] += execution['wall_seconds']
                        costs['GFN2_measured_CPU_seconds'] += execution['cpu_user_seconds'] + execution['cpu_system_seconds']
                result['wall_seconds'] = time.monotonic() - began
                result['artifacts'] = {p.name: pin(p) for p in target.iterdir() if p.is_file()}
                write(target / 'RESULT.json', result)
                return result

            before = costs['GFN2_calls_attempted']
            try:
                with ThreadPoolExecutor(max_workers=2) as workers:
                    lows = list(workers.map(low, ('vacuum', 'ALPB(Water)')))
            finally:
                row['new_GFN2_calls'] = costs['GFN2_calls_attempted'] - before
            row['GFN2'] = {r['task']['medium']: pin(Path(r['task']['directory']) / 'RESULT.json') for r in lows}
            if any(r.get('fatal_process_state') for r in lows):
                stop = 'owned backend process group remained alive'
            require(all(r['status'] == 'accepted' for r in lows), 'GFN2 component unavailable; no retry')
            row.update(status='complete', components={'MACE_eV': mace['energy_eV'],
                       'GFN2_vacuum_hartree': lows[0]['energy_hartree'], 'GFN2_water_hartree': lows[1]['energy_hartree']})
        except Exception as exc:
            row.update(reason=repr(exc), traceback=traceback.format_exc())
        row['wall_seconds'] = time.monotonic() - tick
        write(directory / 'CELL.json', row); row['cell'] = pin(directory / 'CELL.json'); rows.append(row)
        print(json.dumps({'task_id': task['task_id'], 'status': row['status']}), flush=True)
    if native is not None:
        receipt['native_runtime'] = {'versions': native.versions, 'device': native.device,
                                    'actual_model_forwards': native.calls,
                                    'peak_cuda_allocated_bytes': native.torch.cuda.max_memory_allocated()}
        require(native.calls == costs['MACE_calls_attempted'], 'native call counter differs')
    require(costs['MACE_calls_attempted'] <= len(tasks) and costs['GFN2_calls_attempted'] <= 2 * len(tasks),
            'declared component call ceiling exceeded')
    return rows


def summarize(rows):
    cells = {r['task']['task_id']: r for r in rows}; transfers = []; menus = {}
    for source in SOURCES:
        for angle in (0, 90):
            ids = {m: [f'{source}_{shell}_plane{angle}_{m}' for shell in ('inner', 'outer')] for m in METALS}
            missing = [ident for group in ids.values() for ident in group if cells[ident]['status'] != 'complete']
            item = {'source_id': source, 'plane_angle_degree': angle, 'missing_cells': missing,
                    'status': 'unavailable' if missing else 'available', 'own_metal_inner_minus_outer': {}}
            for metal in METALS:
                a, b = ids[metal]
                item['own_metal_inner_minus_outer'][metal] = component_difference(cells[a]['components'], cells[b]['components']) \
                    if all(cells[i]['status'] == 'complete' for i in (a, b)) else None
            if not missing:
                own = item['own_metal_inner_minus_outer']
                item['H_Dy_minus_La_transfer'] = {k: math.fsum([own['Dy'][k], -own['La'][k]]) for k in own['La']}
                h = item['H_Dy_minus_La_transfer']['composite_kcal_mol']
                item['interpretation_at_1_kcal_mol'] = 'near_zero' if abs(h) <= 1. else ('positive' if h > 0 else 'negative')
            transfers.append(item)
        ids = [f'{source}_{pose}_{metal}' for pose in POSES for metal in METALS]
        missing = [i for i in ids if cells[i]['status'] != 'complete']
        menu = {'status': 'unavailable' if missing else 'available', 'candidate_order': list(POSES), 'missing_cells': missing}
        if not missing:
            selected = {}; relative = {}; ties = {}
            for metal in METALS:
                reference = cells[f'{source}_{POSES[0]}_{metal}']['components']
                relative[metal] = {pose: component_difference(cells[f'{source}_{pose}_{metal}']['components'], reference)
                                   for pose in POSES}
                minimum = min(relative[metal][pose]['composite_kcal_mol'] for pose in POSES)
                ties[metal] = [pose for pose in POSES if relative[metal][pose]['composite_kcal_mol'] == minimum]
                selected[metal] = ties[metal][0]
            menu.update(selected=selected, exact_ties=ties, same_metal_relative_energies=relative,
                        selected_Dy_minus_La=component_difference(cells[f"{source}_{selected['Dy']}_Dy"]['components'],
                                                               cells[f"{source}_{selected['La']}_La"]['components']),
                        interpretation='conditional same-water-count energy menu; elemental offset is not affinity')
        menus[source] = menu
    signs = [t['interpretation_at_1_kcal_mol'] for t in transfers if t['status'] == 'available']
    return {'transfers': transfers, 'source_local_common_four_pose_menus': menus,
            'interpretive_resolution_kcal_mol': 1., 'all_four_transfer_contrasts_available': len(signs) == 4,
            'opposite_above_resolution_signs_present': 'positive' in signs and 'negative' in signs,
            'near_zero_present': 'near_zero' in signs, 'no_source_vote_or_population': True,
            'finite_panel_closed_after_declared_remaining_stage': True}


def markdown(result):
    lines = ['# Fixed-source one-water panel', '', f"Stage: {result['stage']}; status: {result['status']}.", '',
             '| Source | Pose | Metal | Availability |', '|---|---|---|---|']
    for row in result['rows']:
        task = row['task']; lines.append(f"| {task['source_id']} | {task['pose_id']} | {task['metal']} | {row['status']} |")
    lines += ['', '| Source | Orientation | La inner−outer | Dy inner−outer | H |', '|---|---:|---:|---:|---:|']
    for transfer in result['summary']['transfers']:
        own = transfer['own_metal_inner_minus_outer']
        vals = [f"{own[m]['composite_kcal_mol']:.6f}" if own[m] is not None else 'unavailable' for m in METALS]
        h = f"{transfer['H_Dy_minus_La_transfer']['composite_kcal_mol']:.6f}" if transfer['status'] == 'available' else 'unavailable'
        lines.append(f"| {transfer['source_id']} | {transfer['plane_angle_degree']} | {vals[0]} | {vals[1]} | {h} |")
    lines += ['', 'Transfer values are kcal/mol; interpretation resolution is1kcal/mol.',
              'All comparisons retain one explicit water. Outer-shell water is not a bulk reference.',
              'Missing cells remain missing; no affinity, population or dry/wet raw-energy selection.', '',
              'Actual new-stage cost: ' + json.dumps(result['actual_cost'], sort_keys=True), '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('scout', 'remaining'), required=True)
    parser.add_argument('--config', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) >= 2,
            'Slurm with at least two allocated CPUs required')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    result = {'protocol_id': PROTOCOL, 'stage': args.stage, 'status': 'execution_failure',
              'config': pin(args.config), 'implementation': pin(__file__), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
              'declared_cells': 16, 'scout_gate_passed': False, 'rows': [],
              'actual_cost': {'MACE_calls_attempted': 0, 'GFN2_calls_attempted': 0, 'model_load_seconds': 0.,
                              'MACE_call_seconds': 0., 'GFN2_solver_wall_seconds': 0., 'GFN2_measured_CPU_seconds': 0.},
              'gradient_calls': 0, 'optimization_steps': 0, 'old256_energy_reuse': False}
    try:
        world = inputs(args.config, args.stage); cfg = world['cfg']
        result.update({k: cfg[k] for k in DESIGN})
        result.update(model=world['mcfg']['model'], GFN2_runtime=world['scfg']['runtime'],
                      GFN2_settings=world['manifest']['settings'], GFN2_tolerances=world['scfg']['tolerances'],
                      native_helper=world['native_helper'], private_atom_count_adaptations=world['adaptations'],
                      preflight=world['preflight'])
        write(out / 'PREFLIGHT.json', {'status': 'passed', 'rows': world['preflight'],
                                     'private_atom_count_adaptations': world['adaptations'], 'energy_calls': 0})
        result['preflight_result'] = pin(out / 'PREFLIGHT.json')
        rows = {t['task_id']: {'task': t, 'status': 'not_yet_scored' if t['admitted'] else 'geometry_rejected'}
                for t in world['tasks']}
        if args.stage == 'remaining':
            prior, reused = previous_scout(world); rows.update(reused)
            result.update(scout_result=cfg['scout_result'], reused_scout_task_ids=list(SCOUT),
                          scout_actual_cost=prior['actual_cost'])
            selected = [t for t in world['tasks'] if t['task_id'] not in SCOUT]
        else:
            selected = [t for t in world['tasks'] if t['task_id'] in SCOUT]
        if args.stage == 'scout' and not all(t['admitted'] for t in selected):
            result.update(status='scout_geometry_unavailable', reason='Declared scout pose rejected; root review required, no alternate')
        else:
            write(out / 'START.json', result)
            executed = execute(world, selected, out, result)
            rows.update({r['task']['task_id']: r for r in executed})
            result['status'] = 'complete' if all(r['status'] in ('complete', 'geometry_rejected') for r in executed) else 'partial'
        result['rows'] = [rows[t['task_id']] for t in world['tasks']]
        result['scout_gate_passed'] = all(rows[i]['status'] == 'complete' for i in SCOUT)
        result['accepted_cells'] = sum(r['status'] == 'complete' for r in result['rows'])
        result['geometry_rejected_cells'] = sum(r['status'] == 'geometry_rejected' for r in result['rows'])
        result['all16_available'] = result['accepted_cells'] == 16
        result['summary'] = summarize(result['rows'])
        result['summary']['finite_panel_closed_after_declared_remaining_stage'] = args.stage == 'remaining'
        result['panel_component_call_totals'] = {key: result['actual_cost'][key] + result.get('scout_actual_cost', {}).get(key, 0)
                                               for key in ('MACE_calls_attempted', 'GFN2_calls_attempted')}
        require(result['panel_component_call_totals']['MACE_calls_attempted'] <= 16 and
                result['panel_component_call_totals']['GFN2_calls_attempted'] <= 32, 'whole finite panel call ceiling exceeded')
        (out / 'REPORT.md').write_text(markdown(result))
        result['report'] = pin(out / 'REPORT.md')
    except Exception as exc:
        result.update(status='execution_failure', reason=repr(exc), traceback=traceback.format_exc())
    finally:
        result.update(wall_seconds=time.monotonic() - started,
                      peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(out / 'RESULT.json', result)
    if result['status'] not in ('complete', 'scout_geometry_unavailable'):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
