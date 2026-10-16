#!/usr/bin/env python3
"""Eight fixed whole-LBT3 standalone GFN2 energy/analytic-gradient calls; Slurm only.

CLI: python LBT3_SOLVENT.py --config CONFIG.json --output FRESH_DIRECTORY
Config pins: implementation, plan, preparation, mace_result, qualification,
backend_helper, gradient_helper, common_helper, python; plus repository,
protocol_id, workers=4, pins=[], runtime and tolerances. Runtime is the exact
{executable, parameters, runtime_libraries, environment} from an accepted
standalone qualification task. Tolerances: charge_closure_e=1e-5,
xtb_json_stdout_hartree=1e-7, engrad_energy_hartree=1e-8,
engrad_coordinates_A=1e-6, xtb_json_population_e=1e-7. Optional mace_forces
pins actual four-cell force result12500. Settings come from qualification12488.
No native ORCA, MACE invocation, geometry response or optimization.
"""
from concurrent.futures import ThreadPoolExecutor
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import socket
import sys
import time
import traceback

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
EXPECTED_SETTINGS = {'accuracy': 0.0001, 'maxiter': 500,
                     'electronic_temperature_K': 300, 'alpb_reference': 'gsolv'}
TOLERANCES = {'charge_closure_e': 1e-5, 'xtb_json_stdout_hartree': 1e-7,
              'engrad_energy_hartree': 1e-8, 'engrad_coordinates_A': 1e-6,
              'xtb_json_population_e': 1e-7}
Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'La': 57, 'Dy': 66}
EV_TO_KCAL = 23.06054783061903  # existing native MACE comparison convention


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return {'path': str(path), 'sha256': digest.hexdigest()}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: ' + record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def xyz(record):
    lines = verify(record).read_text().splitlines()
    rows = [line.split() for line in lines[2:] if line.strip()]
    require(int(lines[0]) == len(rows) == 256 and all(len(row) == 4 for row in rows), 'Whole XYZ atom count/layout')
    positions = [[float(x) for x in row[1:]] for row in rows]
    require(all(math.isfinite(x) for row in positions for x in row), 'Nonfinite XYZ')
    return [row[0] for row in rows], positions


def admit(directory, task, execution, backend, gradients, common, np):
    require(execution['returncode'] == 0 and execution['owned_process_group_empty'], 'Execution failure')
    data = backend.parse_backend(directory, 'xtb', task)
    require(data['status'] == 'available', 'SCC/population failure: ' + str(data.get('failure_reason')))
    energy, charges = data['energy_hartree'], data['charges_e']
    payload = json.loads((directory / 'xtbout.json').read_text())
    require(payload['method'] == 'GFN2-xTB' and payload['xtb version'].startswith('6.7.1'), 'Actual method/version differs')
    require(payload['number of electrons'] == task['valence_electrons']
            and payload['number of unpaired electrons'] == 0, 'Actual effective electronic state differs')
    require(abs(payload['total energy'] - energy) <= TOLERANCES['xtb_json_stdout_hartree'], 'JSON/stdout energy differs')
    require(len(payload['partial charges']) == 256
            and all(math.isfinite(q) for q in payload['partial charges']), 'JSON charge vector invalid')
    require(max(abs(a-b) for a, b in zip(payload['partial charges'], charges))
            <= TOLERANCES['xtb_json_population_e'], 'JSON/charge-file populations differ')
    text = (directory / 'stdout.txt').read_text()
    combined = text + '\n' + (directory / 'stderr.txt').read_text()
    require('normal termination of xtb' in combined.lower(), 'Normal termination absent')
    require(not re.search(r'failed to converge|NOT CONVERGED|abnormal termination|numerical gradient', combined, re.I),
            'SCC failure or unsupported numerical gradient')
    for name, expected in ((r'net charge', -2), (r'unpaired electrons', 0),
                           (r'# electrons', task['valence_electrons']),
                           (r'electronic temp\.', 300), (r':  accuracy', 0.0001)):
        values = re.findall(name + r'\s+([-+\d.Ee]+)', text)
        require(bool(values) and float(values[-1]) == expected, 'Actual standalone setting differs: ' + name)
    for pattern in (r'max\. iterations\s+500', r'restarted\?\s+false', r'PC potential\s+false',
                    r'convergence criteria satisfied after'):
        require(bool(re.search(pattern, text)), 'Actual standalone setting absent: ' + pattern)
    if task['medium'] == 'ALPB(Water)':
        for term in ('Solvation model:               ALPB', 'Reference state                gsolv',
                     'Solvent                        water', 'Grid points                               230',
                     'Ion screening                  false'):
            require(term in text, 'Actual solvent setting absent: ' + term)
    else:
        require('Solvation model:               ALPB' not in text, 'Unexpected solvent')
    require(math.isfinite(energy) and len(charges) == 256 and all(math.isfinite(q) for q in charges), 'Nonfinite energy/populations')
    require(abs(math.fsum(charges) + 2) <= TOLERANCES['charge_closure_e'], 'Charge closure failed')
    candidates = [directory / name for name in ('core.engrad', 'xtb-orca.engrad') if (directory / name).is_file()]
    require(len(candidates) == 1, 'Expected one supported analytic engrad output')
    gradient_pin = pin(candidates[0])
    header = verify(gradient_pin).read_text()
    require('Eh/bohr' in header and 'coordinates in Bohr' in header, 'Analytic gradient units absent')
    grad = gradients.read_engrad(verify(gradient_pin))
    symbols, positions = xyz(task['xyz'])
    require(grad['atom_count'] == 256 and abs(grad['energy_Ha'] - energy) <= TOLERANCES['engrad_energy_hartree'],
            'Gradient count/energy differs')
    require(np.array_equal(grad['atomic_numbers'], [Z[s] for s in symbols]), 'Gradient atom order differs')
    require(np.allclose(grad['coordinates_bohr'] * common.BOHR_TO_A, positions,
                        atol=TOLERANCES['engrad_coordinates_A'], rtol=0), 'Gradient coordinates differ')
    verify(gradient_pin)
    return {'energy_hartree': energy, 'charges_e': charges, 'charge_closure_e': abs(math.fsum(charges) + 2),
            'solver_observation': data, 'xtbout': pin(directory / 'xtbout.json'),
            'gradient': {'artifact': gradient_pin, 'quantity': 'gradient_not_force',
                         'driver': 'standalone xTB 6.7.1 --grad analytic ORCA-format export',
                         'hartree_per_bohr': grad['gradient_Ha_per_bohr'].tolist(),
                         'kcal_mol_A': (grad['gradient_Ha_per_bohr'] * common.HA_TO_KCAL / common.BOHR_TO_A).tolist(),
                         'atom_order_xyz': task['xyz'], 'independent_Dy_force_reference': None}}


def run(config_path, out):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    cfg = json.loads(config_path.read_text())
    require(Path(cfg['repository']).resolve() == ROOT and cfg['protocol_id'], 'Repository/protocol differs')
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Runner differs')
    require(pin(sys.executable)['sha256'] == cfg['python']['sha256'], 'Python differs')
    require(cfg['tolerances'] == TOLERANCES, 'Declared tolerances differ')
    verify(cfg['plan'])
    for record in cfg['pins']:
        verify(record)
    sys.path.insert(0, str(ROOT / 'scripts'))
    import numpy as np
    common = module('lbt3_solvent_common', verify(cfg['common_helper']))
    gradients = module('lbt3_solvent_gradients', verify(cfg['gradient_helper']))
    backend = module('lbt3_solvent_backend', verify(cfg['backend_helper']))
    require(common.BOHR_TO_A == 0.529177210903 and common.HA_TO_KCAL == 627.509474, 'Unit convention changed')
    prep, mace, qualification = [read(cfg[key]) for key in ('preparation', 'mace_result', 'qualification')]
    require(prep['status'] == 'complete' and prep['slurm_job_id'] == '12469', 'Preparation unavailable/different')
    require(mace['status'] == 'complete' and mace['slurm_job_id'] == '12471', 'MACE precedent unavailable/different')
    force_by = {}
    if cfg.get('mace_forces'):
        force_result = read(cfg['mace_forces'])
        require(force_result['status'] == 'complete' and force_result['slurm_job_id'] == '12500'
                and force_result['declared_cells'] == 4 and force_result['analytic_forces_requested']
                and force_result['energy_reference'] == cfg['mace_result']
                and force_result['model'] == mace['model'] and force_result['versions'] == mace['versions'],
                'MACE force precedent/model differs')
        force_by = {row['task_id']: row for row in force_result['rows']}
        require(len(force_by) == len(force_result['rows']) == 4, 'Four unique MACE force rows required')
    require(qualification['status'] == 'passed' and qualification['slurm_job_id'] == '12488'
            and qualification['accepted'] == qualification['attempted'] == 8
            and len(qualification['parity']) == 4 and len(qualification['solvent_component_parity']) == 2
            and all(item['passed'] for item in qualification['parity'] + qualification['solvent_component_parity']),
            'Qualified eight-cell result unavailable/different')
    qcfg, qmanifest = read(qualification['config']), read(qualification['manifest'])
    require(qcfg['staging']['standalone_settings'] == EXPECTED_SETTINGS, 'Qualified settings differ')
    require(cfg['backend_helper'] == qcfg['backend_helper'], 'Execution helper differs from qualification')
    profile = read(qmanifest['profile'])
    runtime = cfg['runtime']
    accepted = [row for row in qualification['rows'] if row['task']['backend'] == 'standalone']
    require(len(accepted) == 4 and all(row['status'] == 'accepted' for row in accepted), 'Standalone qualification incomplete')
    require(all({key: row['task'][key] for key in runtime} == runtime for row in accepted)
            and set(runtime) == {'executable', 'parameters', 'runtime_libraries', 'environment'}, 'Runtime differs from qualification')
    # DY_RUN set this native environment once before launching every subprocess;
    # standalone then overrode its own keys. Reproduce that declared inherited
    # thread/locale/PATH contract without running or importing native ORCA.
    native_environments = [row['task']['environment'] for row in qualification['rows'] if row['task']['backend'] == 'native']
    require(len(native_environments) == 4 and all(env == native_environments[0] for env in native_environments),
            'Qualification inherited environment differs')
    effective_environment = dict(native_environments[0], **runtime['environment'])
    effective_environment.update(OMP_PLACES=None, OMP_PROC_BIND=None)
    runtime_pins = [runtime['executable'], runtime['parameters'], *runtime['runtime_libraries']]
    for record in runtime_pins:
        verify(record)
    bases = [task for row in prep['rows'] for task in row['tasks'] if task['representation'] == 'whole']
    require(len(bases) == 4 and {(t['source_id'], t['metal']) for t in bases}
            == {(s, m) for s in ('7CCO', '7CCN') for m in ('La', 'Dy')}, 'Whole source denominator differs')
    mace_by = {row['task_id']: row for row in mace['rows']}
    signatures = []
    for source in prep['rows']:
        prepared = read(source['preparation'])
        require(source['status'] == 'prepared' and prepared['tasks'] == source['tasks']
                and prepared['sequence'] == 'FIDTNNDGWIEGDELLA'
                and prepared['charge'] == -2 and prepared['protein_charge'] == -5
                and prepared['waters'] == prepared['spectators'] == prepared['caps'] == [], 'Prepared source state differs')
        physical, native, bonds = [read(prepared[key]) for key in ('physical_atoms', 'native_peptide_atoms', 'covalent_bonds')]
        require(len(physical) == 256 and len(native) == 255, 'Prepared atom map differs')
        signatures.append({'atoms': [(a['id'], a['element'], a.get('resname'), a.get('kind')) for a in physical],
                           'charges': [(a['id'], a['charge_e']) for a in native],
                           'bonds': sorted((tuple(sorted((b['atom_a_id'], b['atom_b_id']))), b['forcefield_equilibrium_A']) for b in bonds),
                           'terminal_chemistry': prepared['terminal_chemistry'], 'forcefield': prepared['forcefield'],
                           'pH': prepared['hydrogen_preparation']['pH']})
    require(len(signatures) == 2 and signatures[0] == signatures[1], 'Two-source chemical inventory incompatible')
    source_data = {}
    for task in bases:
        require(task['atom_count'] == 256 and task['metal_index'] == 255 and task['charge'] == -2
                and task['multiplicity'] == (1 if task['metal'] == 'La' else 6)
                and task['source_cap_map'] is None and task['oxidation_state'] == 3, 'Whole source state differs')
        physical = read(task['physical_source'])
        source_data[task['task_id']] = xyz(task['xyz'])
        symbols, positions = source_data[task['task_id']]
        require(symbols == [(task['metal'] if a['element'] == 'M' else a['element']) for a in physical]
                and positions == [a['xyz_A'] for a in physical], 'Whole XYZ/physical map differs')
        require(symbols[255] == task['metal'] and symbols.count(task['metal']) == 1
                and set(symbols[:255]) <= {'H', 'C', 'N', 'O'}, 'Whole source inventory differs')
        old = mace_by[task['task_id']]
        require(old['status'] == 'complete' and all(old[key] == value for key, value in task.items())
                and math.isfinite(old['energy_eV']), 'Saved MACE/source contract differs')
        if force_by:
            force = force_by[task['task_id']]
            require(force['status'] == 'complete' and all(force[key] == value for key, value in task.items())
                    and force['force_convention'] == 'negative derivative of the same native vacuum energy'
                    and force['input_state_check'] == old['input_state_check']
                    and abs(force['energy_eV'] - old['energy_eV']) <= 1e-7, 'Saved MACE force state/energy differs')
            vector = np.asarray(force['forces_eV_A'])
            require(vector.shape == (256, 3) and np.isfinite(vector).all(), 'MACE forces invalid')
    for source in ('7CCO', '7CCN'):
        a, b = [source_data[source + '_whole_' + metal] for metal in ('La', 'Dy')]
        require(a[1] == b[1] and a[0][:255] == b[0][:255], 'Reciprocal source coordinates differ')
    cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
    require(cfg['workers'] == 4 and cpus >= 4, 'Four allocated CPU workers required')
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    receipt = {'schema': 'lbt3_whole_solvent_gradient.v1', 'protocol_id': cfg['protocol_id'],
               'config': pin(config_path), 'plan': cfg['plan'], 'implementation': pin(__file__),
               'preparation': cfg['preparation'], 'mace_result': cfg['mace_result'], 'qualification': cfg['qualification'],
               'mace_forces': cfg.get('mace_forces'),
               'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
               'allocated_cpus': cpus, 'workers': 4, 'declared_endpoints': 8,
               'new_MACE_calls': 0, 'native_ORCA_calls': 0, 'optimizer_steps': 0,
               'status': 'running', 'rows': [], 'corrections': [],
               'gradient_units': {'hartree_to_kcal_mol': common.HA_TO_KCAL, 'bohr_to_A': common.BOHR_TO_A,
                                  'eV_to_kcal_mol': EV_TO_KCAL, 'hartree_to_eV': common.HA_TO_KCAL / EV_TO_KCAL},
               'scope': 'fixed whole-peptide composite; no affinity, Mex composite, or independent force qualification'}
    write(out / 'START.json', receipt)
    tasks = []
    try:
        for base in bases:
            symbols, _ = source_data[base['task_id']]
            electrons = math.fsum(math.fsum(profile['element'][s]['refocc']) for s in symbols) + 2
            require(electrons == int(electrons) and int(electrons) % 2 == 0, 'Invalid GFN2 effective electron count')
            for medium in ('vacuum', 'ALPB(Water)'):
                task_id = base['task_id'] + '__' + ('water' if medium != 'vacuum' else 'vacuum')
                directory = out / task_id
                directory.mkdir()
                with (directory / 'core.xyz').open('xb') as stream:
                    stream.write(verify(base['xyz']).read_bytes())
                command = [runtime['executable']['path'], 'core.xyz', '--gfn', '2', '--chrg', '-2', '--uhf', '0',
                           '--etemp', '300', '--acc', '0.0001', '--iterations', '500', '--norestart', '--json', '--pop', '--grad']
                if medium != 'vacuum':
                    command += ['--alpb', 'water', 'gsolv']
                task = dict(base, task_id=task_id, origin_task_id=base['task_id'], source_xyz=base['xyz'],
                            xyz=pin(directory / 'core.xyz'), directory=str(directory), medium=medium,
                            physical_multiplicity=base['multiplicity'], multiplicity=1, uhf=0,
                            physical_electronic_state=base['electronic_state'],
                            electronic_state='original GFN2 spin-independent f-in-core effective singlet', spin_multiplicity=1,
                            valence_electrons=int(electrons), command=command,
                            state_convention='GFN2 original spin-independent f-in-core; physical spin retained separately')
                write(directory / 'TASK.json', task)
                tasks.append(dict(task, task_record=pin(directory / 'TASK.json')))
        write(out / 'MANIFEST.json', {'tasks': tasks, 'runtime': runtime, 'effective_environment': effective_environment, 'settings': EXPECTED_SETTINGS,
                                      'tolerances': TOLERANCES, 'config': pin(config_path)})
        receipt['manifest'] = pin(out / 'MANIFEST.json')

        def execute(task):
            directory = Path(task['directory'])
            row = {'task': task, 'status': 'unavailable', 'attempted': False,
                   'energy_hartree': None, 'charges_e': None, 'gradient': None}
            try:
                for record in [task['task_record'], task['xyz'], *runtime_pins]:
                    verify(record)
                env = os.environ.copy()
                for key, value in effective_environment.items():
                    if value is None:
                        env.pop(key, None)
                    else:
                        env[key] = value
                for key in ('OMP_PLACES', 'OMP_PROC_BIND'):
                    env.pop(key, None)
                row['attempted'] = True
                execution = backend.call(task['command'], directory, env)
                row['execution'] = pin(directory / 'EXECUTION.json')
                verify(task['xyz'])
                row.update(admit(directory, task, execution, backend, gradients, common, np), status='accepted')
                for record in runtime_pins:
                    verify(record)
            except Exception as exc:
                row.update(status='unavailable', failure=repr(exc), traceback=traceback.format_exc())
            row['artifacts'] = {path.name: pin(path) for path in directory.iterdir() if path.is_file()}
            write(directory / 'RESULT.json', row)
            print(json.dumps({'id': task['task_id'], 'status': row['status'], 'failure': row.get('failure')}), flush=True)
            return row

        with ThreadPoolExecutor(max_workers=4) as pool:
            rows = list(pool.map(execute, tasks))
        receipt['rows'] = rows
        by = {(row['task']['origin_task_id'], row['task']['medium']): row for row in rows}
        for base in bases:
            a, b = [by[(base['task_id'], medium)] for medium in ('ALPB(Water)', 'vacuum')]
            item = {'origin_task_id': base['task_id'], 'source_id': base['source_id'], 'metal': base['metal'],
                    'xyz': base['xyz'], 'status': 'unavailable', 'solvent_transfer_hartree': None,
                    'solvent_gradient_kcal_mol_A': None, 'composite_gradient': None,
                    'MACE_forces_available': bool(force_by), 'composite_forces_eV_A': None, 'composite_kcal_mol': None,
                    'MACE_vacuum_eV': mace_by[base['task_id']]['energy_eV']}
            if a['status'] == b['status'] == 'accepted':
                delta = math.fsum([a['energy_hartree'], -b['energy_hartree']])
                dg = np.asarray(a['gradient']['hartree_per_bohr']) - np.asarray(b['gradient']['hartree_per_bohr'])
                item.update(status='available', solvent_transfer_hartree=delta,
                            solvent_transfer_kcal_mol=delta * common.HA_TO_KCAL,
                            solvent_gradient_hartree_per_bohr=dg.tolist(),
                            solvent_gradient_kcal_mol_A=(dg * common.HA_TO_KCAL / common.BOHR_TO_A).tolist(),
                            gradient_quantity='ALPB-minus-vacuum gradient; not force or full composite gradient',
                            composite_kcal_mol=math.fsum([item['MACE_vacuum_eV'] * EV_TO_KCAL, delta * common.HA_TO_KCAL]))
                if force_by:
                    mace_force = np.asarray(force_by[base['task_id']]['forces_eV_A'])
                    composite_force = mace_force - dg * (common.HA_TO_KCAL / EV_TO_KCAL) / common.BOHR_TO_A
                    require(np.isfinite(composite_force).all(), 'Composite force nonfinite')
                    item.update(MACE_forces_eV_A=mace_force.tolist(), composite_forces_eV_A=composite_force.tolist(),
                                composite_gradient={'kcal_mol_A': (-composite_force * EV_TO_KCAL).tolist(),
                                                    'quantity': 'gradient of MACE_vacuum plus GFN2_water minus GFN2_vacuum'},
                                composite_force_formula='F_MACE - (g_GFN2_water - g_GFN2_vacuum) * hartree_to_eV / bohr_to_A',
                                independent_composite_directional_derivative_check=False)
            receipt['corrections'].append(item)
        corrections = {(item['source_id'], item['metal']): item for item in receipt['corrections']}
        contrasts = []
        for source in ('7CCO', '7CCN'):
            la, dy = [corrections[(source, metal)] for metal in ('La', 'Dy')]
            vacuum = math.fsum([dy['MACE_vacuum_eV'], -la['MACE_vacuum_eV']]) * EV_TO_KCAL
            item = {'source_id': source, 'sign': 'Dy_minus_La', 'vacuum_kcal_mol': vacuum,
                    'solvent_contribution_kcal_mol': None, 'composite_kcal_mol': None, 'status': 'unavailable'}
            if la['status'] == dy['status'] == 'available':
                delta = math.fsum([dy['solvent_transfer_hartree'], -la['solvent_transfer_hartree']]) * common.HA_TO_KCAL
                item.update(status='available', solvent_contribution_kcal_mol=delta,
                            composite_kcal_mol=math.fsum([vacuum, delta]))
            contrasts.append(item)
        receipt['source_contrasts'] = contrasts
        receipt['source_contrast_difference'] = {'sign': '7CCN_minus_7CCO', 'vacuum_kcal_mol': math.fsum([contrasts[1]['vacuum_kcal_mol'], -contrasts[0]['vacuum_kcal_mol']]),
                                                 'solvent_kcal_mol': None, 'composite_kcal_mol': None}
        if all(item['status'] == 'available' for item in contrasts):
            for key, field in (('solvent_kcal_mol', 'solvent_contribution_kcal_mol'), ('composite_kcal_mol', 'composite_kcal_mol')):
                receipt['source_contrast_difference'][key] = math.fsum([contrasts[1][field], -contrasts[0][field]])
        pools = {'interpretation': 'conditional fixed-inventory two-source minima; no population or affinity',
                 'source_order': ['7CCO', '7CCN'], 'metal_order': ['La', 'Dy'], 'chemical_inventory_matches': True}
        for name in ('vacuum', 'composite'):
            matrix = [[corrections[(source, metal)]['MACE_vacuum_eV'] * EV_TO_KCAL if name == 'vacuum'
                       else corrections[(source, metal)]['composite_kcal_mol'] for source in ('7CCO', '7CCN')] for metal in ('La', 'Dy')]
            entry = {'status': 'unavailable', 'matrix_kcal_mol': matrix, 'selected_sources': None, 'Dy_minus_La_kcal_mol': None}
            if all(value is not None for row in matrix for value in row):
                minima = [min(row) for row in matrix]
                selected = {metal: [source for source, value in zip(('7CCO', '7CCN'), row) if value == low]
                            for metal, row, low in zip(('La', 'Dy'), matrix, minima)}
                entry.update(status='available', minimum_energy_kcal_mol=minima, selected_sources=selected,
                             Dy_minus_La_kcal_mol=math.fsum([minima[1], -minima[0]]))
            pools[name] = entry
        if pools['composite']['status'] == 'available':
            pools['selection_changed'] = {metal: pools['vacuum']['selected_sources'][metal] != pools['composite']['selected_sources'][metal] for metal in ('La', 'Dy')}
            pools['contrast_change_kcal_mol'] = math.fsum([pools['composite']['Dy_minus_La_kcal_mol'], -pools['vacuum']['Dy_minus_La_kcal_mol']])
        receipt['common_two_source_pool'] = pools
        receipt.update(status='complete' if all(row['status'] == 'accepted' for row in rows) else 'partial',
                       attempted=sum(row['attempted'] for row in rows), accepted=sum(row['status'] == 'accepted' for row in rows))
        verify(receipt['config'])
    except BaseException as exc:
        receipt.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    receipt['elapsed_seconds'] = time.monotonic() - started
    receipt['application_allocated_CPU_seconds'] = receipt['elapsed_seconds'] * cpus
    receipt['allocation'] = {key: os.environ.get(key) for key in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_CPU', 'SLURM_MEM_PER_NODE')}
    write(out / 'RESULT.json', receipt)
    return receipt['status'] == 'complete'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not run(args.config.resolve(), args.output.resolve()):
        raise SystemExit(1)
