#!/usr/bin/env python3
"""Pinned LBT3 finite-direction energies and saved derivative comparison.

CLI: LBT3_FD.py {gpu,cpu,collect} --config CONFIG --output FRESH_DIRECTORY
Common config: protocol_id, repository, implementation/plan/profile_result/
solvent_result pins, pins[], workers=4, tolerances={absolute_kcal_mol_A:0.2,
relative:0.05}. collect additionally requires stage_results={gpu:pin,cpu:pin}.
Model, helpers, runtime, settings and Python are inherited from accepted12510
and12500. Exactly eight MACE and sixteen standalone GFN2 fixed energies; no
gradients requested by either FD execution stage, no geometry generation.
"""
from concurrent.futures import ThreadPoolExecutor
import argparse
import hashlib
import importlib.metadata
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
EV_TO_KCAL = 23.06054783061903
HA_TO_KCAL = 627.509474
STEPS = (0.001, 0.0005)
METALS = ('La', 'Dy')
STATES = ('La', 'Dy', 'common', 'Dy_minus_La')
COMPONENTS = ('MACE_vacuum', 'GFN2_vacuum', 'GFN2_water', 'solvent', 'composite')
TOLERANCES = {'absolute_kcal_mol_A': 0.2, 'relative': 0.05}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve(); digest = hashlib.sha256()
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
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def xyz(record):
    lines = verify(record).read_text().splitlines(); rows = [line.split() for line in lines[2:] if line.strip()]
    require(int(lines[0]) == len(rows) == 256 and all(len(row) == 4 for row in rows), 'Whole XYZ layout differs')
    symbols = [row[0] for row in rows]; positions = [[float(x) for x in row[1:]] for row in rows]
    require(all(math.isfinite(x) for row in positions for x in row), 'Nonfinite coordinates')
    return symbols, positions


class NativeMACE:
    """Existing native warm-model call; no coordinates, optimizer or score policy."""
    def __init__(self, config):
        require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
        require(pin(sys.executable)['sha256'] == config['python']['sha256'], 'Native MACE Python differs')
        for record in config['helper_pins'].values():
            verify(record)
        import torch
        from ase import Atoms
        from ase.data import atomic_numbers
        from mace.calculators import mace_omol
        sys.path.insert(0, str(verify(config['helper_pins']['mace_omol.py']).parent))
        from mace_omol import input_batch
        self.versions = {key: importlib.metadata.version(key) for key in ('torch', 'mace-torch', 'ase')}
        require(self.versions == config['expected_versions'], 'Native MACE software differs')
        require(torch.cuda.is_available() and torch.cuda.device_count() == 1
                and os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'), 'One allocated GPU required')
        torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK'])); torch.set_default_dtype(torch.float64)
        torch.cuda.reset_peak_memory_stats()
        model_config = config['model']
        require(model_config['dtype'] == 'float64' and model_config['execution_adapter'] is None
                and model_config['solvent'] is None and not model_config['extra_dispersion'], 'Native recipe differs')
        self.calc = mace_omol(model=str(verify(model_config['checkpoint'])), device='cuda', default_dtype='float64')
        self.model = self.calc.models[0]
        require(type(self.model).__name__ == 'ScaleShiftMACE'
                and dict(self.model.embedding_specs) == model_config['embedding_specs']
                and list(self.model.heads) == ['omol'] and self.calc.head == 'omol'
                and self.calc.energy_units_to_eV == 1., 'Native checkpoint/head/units differ')
        for parameter in self.model.parameters():
            parameter.requires_grad_(False)
        self.parameter_versions = [parameter._version for parameter in self.model.parameters()]
        self.torch, self.Atoms, self.Z, self.input_batch = torch, Atoms, atomic_numbers, input_batch
        self.calls = 0; self.device = torch.cuda.get_device_name()

    def evaluate(self, task, symbols, positions, compute_force=False):
        tick = time.monotonic(); torch = self.torch
        require(len(symbols) == len(positions) == 256 and symbols[255] == task['metal']
                and symbols.count(task['metal']) == 1 and task['metal_index'] == 255
                and task['charge'] == -2 and task['multiplicity'] == (1 if task['metal'] == 'La' else 6), 'MACE state differs')
        require(all(len(row) == 3 and all(math.isfinite(x) for x in row) for row in positions), 'Invalid MACE coordinates')
        require(set(symbols) <= {'H', 'C', 'N', 'O', task['metal']}
                and {self.Z[s] for s in symbols} <= set(self.model.atomic_numbers.tolist()), 'MACE element unsupported')
        require((sum(self.Z[s] for s in symbols) - task['charge'] - task['multiplicity'] + 1) % 2 == 0,
                'All-electron physical spin parity differs')
        atoms = self.Atoms(symbols, positions=positions, pbc=False)
        atoms.info.update(charge=task['charge'], spin=task['multiplicity'])
        batch = self.calc._atoms_to_batch(atoms)
        state = self.input_batch(self.calc, atoms, task['charge'], task['multiplicity'], batch, 255)
        with torch.enable_grad() if compute_force else torch.no_grad():
            self.calls += 1
            value = self.model(batch.to_dict(), training=False, compute_force=compute_force,
                               compute_virials=False, compute_stress=False, compute_displacement=False,
                               compute_hessian=False, compute_edge_forces=False, compute_atomic_stresses=False)
        energy = float(value['energy'].detach().cpu().item()); forces = None
        require(math.isfinite(energy), 'Nonfinite MACE energy')
        require([p._version for p in self.model.parameters()] == self.parameter_versions, 'Model parameters changed')
        if compute_force:
            force = value['forces'].detach().cpu()
            require(tuple(force.shape) == (256, 3) and bool(torch.isfinite(force).all()), 'Invalid native force')
            forces = force.tolist()
        else:
            require(value['forces'] is None, 'Unexpected force output')
        return {'energy_eV': energy, 'forces_eV_A': forces, 'input_state_check': state,
                'force_convention': 'negative derivative of the same native vacuum energy' if compute_force else None,
                'wall_seconds': time.monotonic() - tick}


def inputs(config_path, stage):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'All FD stages require Slurm')
    cfg = json.loads(config_path.read_text())
    require(Path(cfg['repository']).resolve() == ROOT and cfg['protocol_id'], 'Repository/protocol differs')
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Implementation changed')
    require(cfg['workers'] == 4 and cfg['tolerances'] == TOLERANCES, 'Declared settings changed')
    verify(cfg['plan'])
    for record in cfg['pins']:
        verify(record)
    profile, solvent = read(cfg['profile_result']), read(cfg['solvent_result'])
    require(profile['status'] == 'complete' and profile['solvent_result'] == cfg['solvent_result']
            and profile['plan'] == cfg['plan'], 'Completed profile/plan binding differs')
    require(solvent['status'] == 'complete' and solvent['accepted'] == solvent['attempted'] == 8
            and solvent['slurm_job_id'] == '12510', 'Whole origin solvent unavailable/different')
    finite = read(profile['finite_check']); scfg = read(solvent['config'])
    force_result = read(solvent['mace_forces']); mcfg = read(force_result['config'])
    require(force_result['status'] == 'complete' and force_result['slurm_job_id'] == '12500'
            and force_result['model'] == mcfg['model'] and force_result['versions'] == mcfg['expected_versions'],
            'Native force/model lineage differs')
    require(pin(sys.executable)['sha256'] == (mcfg if stage == 'gpu' else scfg)['python']['sha256'], 'Stage Python differs')
    require(finite['status'] == 'prepared' and finite['source_id'] == '7CCO'
            and finite['source_guard']['pass'] and finite['kinematics_check_at_positive_half_step']['pass'],
            'Finite geometry preparation unavailable')
    require(finite['numerical_admission_from_plan'] == dict(TOLERANCES, steps_A=list(STEPS)), 'Finite tolerance differs')
    analytic = finite['analytic_energy_derivatives_kcal_mol_A']
    require(set(analytic) == set(COMPONENTS) and all(set(analytic[c]) == set(STATES) for c in COMPONENTS)
            and all(math.isfinite(value) for values in analytic.values() for value in values.values()), 'Analytic projections incomplete')
    tasks = finite['tasks']; geometries = {g['geometry_id']: g for g in finite['geometries']}
    require(len(tasks) == len({t['task_id'] for t in tasks}) == 8 and len(geometries) == 4
            and {(t['metal'], t['signed_step_A']) for t in tasks} == {(m, sign*h) for m in METALS for h in STEPS for sign in (-1, 1)},
            'Eight finite states changed')
    origins = {row['task_id']: row for row in force_result['rows']}; coordinates = {}
    for task in tasks:
        require(task['status'] == 'prepared' and task['source_id'] == '7CCO'
                and task['representation'] == 'whole' and task['source_cap_map'] is None, 'Candidate state unavailable')
        geometry = read(task['geometry'])
        require(geometry == geometries[task['geometry_id']] and geometry['status'] == 'prepared'
                and geometry['inherited_guard']['pass'] and geometry['signed_step_A'] == task['signed_step_A']
                and geometry['mode_index'] == finite['mode_index'] and geometry['mode_name'] == finite['mode']['name'],
                'Candidate geometry binding differs')
        origin = origins[task['origin_task_id']]
        require(origin['xyz'] == task['origin_xyz'] and origin['physical_source'] == task['origin_physical_source']
                and all(origin[key] == task[key] for key in ('source_id', 'metal', 'charge', 'multiplicity', 'atom_count', 'metal_index', 'oxidation_state')),
                'Origin physical state changed')
        symbols, positions = xyz(task['xyz']); physical = read(task['physical_source'])
        require(len(physical) == 256 and symbols == [(task['metal'] if a['element'] == 'M' else a['element']) for a in physical]
                and positions == [a['xyz_A'] for a in physical], 'Candidate XYZ/physical map differs')
        coordinates[task['task_id']] = (symbols, positions)
    for step in (-0.001, -0.0005, 0.0005, 0.001):
        pair = [next(t for t in tasks if t['metal'] == metal and t['signed_step_A'] == step) for metal in METALS]
        a, b = [coordinates[t['task_id']] for t in pair]
        require(a[1] == b[1] and a[0][:255] == b[0][:255], 'Paired candidate coordinates differ')
    design = {key: cfg[key] for key in ('plan', 'profile_result', 'solvent_result', 'tolerances')}
    design.update(implementation_sha256=cfg['implementation']['sha256'], finite_check=profile['finite_check'])
    return cfg, design, finite, tasks, coordinates, solvent, scfg, mcfg


def admit_energy(directory, task, execution, backend, tolerances):
    require(execution['returncode'] == 0 and execution['owned_process_group_empty'], 'Standalone execution failed')
    data = backend.parse_backend(directory, 'xtb', task)
    require(data['status'] == 'available', 'SCC/population failure: ' + str(data.get('failure_reason')))
    energy, charges = data['energy_hartree'], data['charges_e']
    payload = json.loads((directory / 'xtbout.json').read_text())
    require(payload['method'] == 'GFN2-xTB' and payload['xtb version'].startswith('6.7.1')
            and payload['number of electrons'] == task['valence_electrons']
            and payload['number of unpaired electrons'] == 0, 'Standalone actual method/state differs')
    require(abs(payload['total energy'] - energy) <= tolerances['xtb_json_stdout_hartree'], 'JSON energy differs')
    require(len(charges) == len(payload['partial charges']) == 256
            and all(math.isfinite(q) for q in charges + payload['partial charges'])
            and abs(math.fsum(charges) + 2) <= tolerances['charge_closure_e']
            and max(abs(a-b) for a, b in zip(charges, payload['partial charges'])) <= tolerances['xtb_json_population_e'], 'Population checks failed')
    text = (directory / 'stdout.txt').read_text(); combined = text + '\n' + (directory / 'stderr.txt').read_text()
    require('normal termination of xtb' in combined.lower()
            and not re.search(r'failed to converge|NOT CONVERGED|abnormal termination', combined, re.I), 'Standalone termination/convergence failed')
    for name, expected in ((r'net charge', -2), (r'unpaired electrons', 0), (r'# electrons', task['valence_electrons']),
                           (r'electronic temp\.', 300), (r':  accuracy', 0.0001)):
        values = re.findall(name + r'\s+([-+\d.Ee]+)', text)
        require(bool(values) and float(values[-1]) == expected, 'Standalone actual setting differs: ' + name)
    for pattern in (r'max\. iterations\s+500', r'restarted\?\s+false', r'PC potential\s+false', r'convergence criteria satisfied after'):
        require(bool(re.search(pattern, text)), 'Standalone actual setting absent: ' + pattern)
    if task['medium'] == 'ALPB(Water)':
        for term in ('Solvation model:               ALPB', 'Reference state                gsolv', 'Solvent                        water',
                     'Grid points                               230', 'Ion screening                  false'):
            require(term in text, 'Standalone solvent setting differs: ' + term)
    else:
        require('Solvation model:               ALPB' not in text, 'Unexpected solvent')
    return dict(energy_hartree=energy, charges_e=charges, solver_observation=data,
                xtbout=pin(directory / 'xtbout.json'), gradient_requested=False)


def gpu_run(tasks, coordinates, mcfg, out, receipt):
    evaluator = NativeMACE(mcfg); stop = None; rows = []
    for task in tasks:
        row = {'task': task, 'status': 'unavailable', 'call_attempted': False, 'energy_eV': None}
        if stop:
            row['failure'] = 'Not attempted after model failure: ' + stop
        else:
            try:
                verify(task['xyz']); row['call_attempted'] = True
                row.update(evaluator.evaluate(task, *coordinates[task['task_id']], compute_force=False), status='accepted')
                verify(task['xyz'])
            except Exception as exc:
                stop = repr(exc); row.update(status='unavailable', failure=stop, traceback=traceback.format_exc())
        rows.append(row); write(out / (task['task_id'] + '.json'), row)
        print(json.dumps({'id': task['task_id'], 'status': row['status']}), flush=True)
    receipt.update(rows=rows, molecular_calls_attempted=evaluator.calls, accepted=sum(r['status'] == 'accepted' for r in rows),
                   versions=evaluator.versions, device=evaluator.device,
                   peak_cuda_allocated_bytes=evaluator.torch.cuda.max_memory_allocated(), model=mcfg['model'])


def cpu_run(tasks, solvent, scfg, out, receipt):
    require(int(os.environ['SLURM_CPUS_PER_TASK']) >= 4, 'Four allocated CPU workers required')
    backend = module('lbt3_fd_backend', verify(scfg['backend_helper']))
    manifest = read(solvent['manifest']); runtime = manifest['runtime']; envspec = manifest['effective_environment']
    require(runtime == scfg['runtime'] and manifest['settings'] == {'accuracy': 0.0001, 'maxiter': 500,
            'electronic_temperature_K': 300, 'alpb_reference': 'gsolv'}, 'Qualified standalone recipe differs')
    runtime_pins = [runtime['executable'], runtime['parameters'], *runtime['runtime_libraries']]
    for record in runtime_pins:
        verify(record)
    origin_rows = {(r['task']['origin_task_id'], r['task']['medium']): r for r in solvent['rows']}
    work = []
    for task in tasks:
        for medium in ('vacuum', 'ALPB(Water)'):
            original = origin_rows[(task['origin_task_id'], medium)]
            require(original['status'] == 'accepted' and original['task']['source_xyz'] == task['origin_xyz'], 'GFN2 origin lineage differs')
            identity = task['task_id'] + '__' + ('water' if medium != 'vacuum' else 'vacuum')
            directory = out / identity; directory.mkdir()
            with (directory / 'core.xyz').open('xb') as stream:
                stream.write(verify(task['xyz']).read_bytes())
            command = [runtime['executable']['path'], 'core.xyz', '--gfn', '2', '--chrg', '-2', '--uhf', '0',
                       '--etemp', '300', '--acc', '0.0001', '--iterations', '500', '--norestart', '--json', '--pop', '--sp']
            if medium != 'vacuum':
                command += ['--alpb', 'water', 'gsolv']
            local = dict(task, task_id=identity, finite_task_id=task['task_id'], physical_multiplicity=task['multiplicity'],
                         multiplicity=1, spin_multiplicity=1, uhf=0, medium=medium, command=command,
                         physical_electronic_state=task['electronic_state'], electronic_state='GFN2 spin-independent f-in-core effective singlet',
                         valence_electrons=original['task']['valence_electrons'], source_xyz=task['xyz'],
                         xyz=pin(directory / 'core.xyz'), directory=str(directory), gradient_requested=False)
            write(directory / 'TASK.json', local); work.append(dict(local, task_record=pin(directory / 'TASK.json')))
    write(out / 'MANIFEST.json', {'tasks': work, 'runtime': runtime, 'effective_environment': envspec,
                                  'settings': manifest['settings'], 'design': receipt['design']})
    receipt['manifest'] = pin(out / 'MANIFEST.json')

    def execute(task):
        directory = Path(task['directory']); row = {'task': task, 'status': 'unavailable', 'call_attempted': False, 'energy_hartree': None}
        try:
            for record in [task['task_record'], task['xyz'], *runtime_pins]:
                verify(record)
            env = os.environ.copy()
            for key, value in envspec.items():
                if value is None:
                    env.pop(key, None)
                else:
                    env[key] = value
            row['call_attempted'] = True
            execution = backend.call(task['command'], directory, env); row['execution'] = pin(directory / 'EXECUTION.json')
            verify(task['xyz'])
            row.update(admit_energy(directory, task, execution, backend, scfg['tolerances']), status='accepted')
            for record in runtime_pins:
                verify(record)
        except Exception as exc:
            row.update(status='unavailable', failure=repr(exc), traceback=traceback.format_exc())
        row['artifacts'] = {path.name: pin(path) for path in directory.iterdir() if path.is_file()}
        write(directory / 'RESULT.json', row)
        print(json.dumps({'id': task['task_id'], 'status': row['status']}), flush=True)
        return row

    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(execute, work))
    receipt.update(rows=rows, molecular_calls_attempted=sum(r['call_attempted'] for r in rows),
                   accepted=sum(r['status'] == 'accepted' for r in rows))


def collect(cfg, design, finite, tasks, receipt):
    require(set(cfg['stage_results']) == {'gpu', 'cpu'}, 'Both stage result pins required')
    stages = {stage: read(record) for stage, record in cfg['stage_results'].items()}
    for stage, result in stages.items():
        require(result['stage'] == stage and result['design'] == design, 'Stage design binding differs')
        require(result['status'] in ('complete', 'partial'), 'Stage execution failure is not admissible')
        stage_config = read(result['config']); verify(result['implementation'])
        require(all(stage_config[key] == cfg[key] for key in ('plan', 'profile_result', 'solvent_result', 'tolerances'))
                and stage_config['implementation']['sha256'] == cfg['implementation']['sha256'], 'Stage config binding differs')
    gpu_rows = {row['task']['task_id']: row for row in stages['gpu']['rows']}
    cpu_rows = {(row['task']['finite_task_id'], row['task']['medium']): row for row in stages['cpu']['rows']}
    require(len(gpu_rows) == len(stages['gpu']['rows']) and len(cpu_rows) == len(stages['cpu']['rows']), 'Duplicate stage rows')
    expected_ids = {task['task_id'] for task in tasks}
    require(set(gpu_rows) <= expected_ids and set(cpu_rows) <= {(ident, medium) for ident in expected_ids
            for medium in ('vacuum', 'ALPB(Water)')}, 'Undeclared stage rows')
    values = {}; availability = []
    for task in tasks:
        ident = task['task_id']; row = gpu_rows.get(ident)
        if row is not None:
            require(row['task'] == task, 'MACE task differs from finite preparation')
        mace = row['energy_eV'] if row is not None and row['status'] == 'accepted' else None
        energies = {'MACE_vacuum': mace}
        for medium, component in (('vacuum', 'GFN2_vacuum'), ('ALPB(Water)', 'GFN2_water')):
            row = cpu_rows.get((ident, medium))
            if row is not None:
                require(row['task']['source_xyz'] == task['xyz'] and row['task']['charge'] == -2
                        and row['task']['physical_multiplicity'] == task['multiplicity'], 'GFN2 candidate state differs')
            energies[component] = row['energy_hartree'] if row is not None and row['status'] == 'accepted' else None
        require(all(value is None or math.isfinite(value) for value in energies.values()), 'Nonfinite admitted endpoint')
        values[(task['metal'], task['signed_step_A'])] = energies
        availability.append({'task_id': ident, 'metal': task['metal'], 'signed_step_A': task['signed_step_A'],
                             'MACE_eV': energies['MACE_vacuum'], 'GFN2_vacuum_hartree': energies['GFN2_vacuum'],
                             'GFN2_water_hartree': energies['GFN2_water'], 'complete': all(v is not None for v in energies.values())})
    derivatives = []; step_values = {}
    for h in STEPS:
        component_values = {component: {} for component in COMPONENTS}
        for metal in METALS:
            plus, minus = values[(metal, h)], values[(metal, -h)]
            for component in COMPONENTS[:3]:
                a, b = plus[component], minus[component]
                factor = EV_TO_KCAL if component == 'MACE_vacuum' else HA_TO_KCAL
                component_values[component][metal] = None if a is None or b is None else math.fsum([a, -b]) * factor / (2*h)
            a, b = component_values['GFN2_water'][metal], component_values['GFN2_vacuum'][metal]
            component_values['solvent'][metal] = None if a is None or b is None else math.fsum([a, -b])
            a, b = component_values['MACE_vacuum'][metal], component_values['solvent'][metal]
            component_values['composite'][metal] = None if a is None or b is None else math.fsum([a, b])
        for component, states in component_values.items():
            la, dy = states['La'], states['Dy']
            states['common'] = None if la is None or dy is None else math.fsum([la, dy]) / 2
            states['Dy_minus_La'] = None if la is None or dy is None else math.fsum([dy, -la])
            for state, numerical in states.items():
                analytic = finite['analytic_energy_derivatives_kcal_mol_A'][component][state]
                threshold = TOLERANCES['absolute_kcal_mol_A'] + TOLERANCES['relative'] * abs(analytic)
                residual = None if numerical is None else numerical - analytic
                derivatives.append({'step_A': h, 'component': component, 'state': state,
                    'analytic_kcal_mol_A': analytic, 'numerical_kcal_mol_A': numerical,
                    'residual_kcal_mol_A': residual, 'tolerance_kcal_mol_A': threshold,
                    'within_declared_resolution': residual is not None and abs(residual) <= threshold,
                    'acceptance_role': 'qualification' if component == 'composite' else 'component diagnostic'})
                step_values[(h, component, state)] = numerical
    step_comparisons = []
    for component in COMPONENTS:
        for state in STATES:
            coarse, fine = [step_values[(h, component, state)] for h in STEPS]
            step_comparisons.append({'component': component, 'state': state,
                'fine_minus_coarse_kcal_mol_A': None if coarse is None or fine is None else math.fsum([fine, -coarse])})
    composite = [row for row in derivatives if row['component'] == 'composite']
    own = [row for row in composite if row['state'] in METALS]
    passed = all(row['within_declared_resolution'] for row in composite)
    complete = all(row['complete'] for row in availability)
    receipt.update(status='passed' if complete and passed else ('failed_derivative' if complete else 'unavailable'),
                   stage_results=cfg['stage_results'], endpoint_availability=availability,
                   derivatives=derivatives, step_comparisons=step_comparisons,
                   individual_metal_composite_passed=all(row['within_declared_resolution'] for row in own),
                   composite_passed_all_states=passed, molecular_calls_attempted=0,
                   originating_call_counts={stage: result.get('molecular_calls_attempted', 0) for stage, result in stages.items()},
                   formula='central (E(+h)-E(-h))/(2h); h is signed maximum atom chord displacement in Å',
                   residual_tolerance='0.2 + 0.05*abs(saved analytic derivative), kcal/mol/Å; both steps retained',
                   limits='One force-nominated direction; no affinity, global force validation, optimizer or biological label')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('gpu', 'cpu', 'collect'))
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); cfg, design, finite, tasks, coords, solvent, scfg, mcfg = inputs(args.config.resolve(), args.stage)
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
    receipt = {'schema': 'lbt3_fd_energy_or_comparison.v1', 'stage': args.stage, 'config': pin(args.config),
               'implementation': pin(__file__), 'design': design, 'slurm_job_id': os.environ['SLURM_JOB_ID'],
               'host': socket.gethostname(), 'status': 'running', 'rows': [], 'molecular_calls_attempted': 0,
               'declared_calls': {'gpu': 8, 'cpu': 16, 'collect': 0}[args.stage], 'new_gradients_requested': 0,
               'optimizer_steps': 0, 'allocated_cpus': cpus}
    write(out / 'START.json', receipt)
    try:
        if args.stage == 'gpu':
            gpu_run(tasks, coords, mcfg, out, receipt)
        elif args.stage == 'cpu':
            cpu_run(tasks, solvent, scfg, out, receipt)
        else:
            collect(cfg, design, finite, tasks, receipt)
        if args.stage != 'collect':
            receipt['status'] = 'complete' if receipt.get('accepted') == receipt['declared_calls'] else 'partial'
        verify(receipt['config'])
    except BaseException as exc:
        receipt.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    receipt['elapsed_seconds'] = time.monotonic() - started
    receipt['application_allocated_CPU_seconds'] = receipt['elapsed_seconds'] * cpus
    receipt['application_allocated_GPU_seconds'] = receipt['elapsed_seconds'] if args.stage == 'gpu' else 0
    receipt['allocation'] = {key: os.environ.get(key) for key in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_CPU', 'SLURM_MEM_PER_NODE', 'CUDA_VISIBLE_DEVICES')}
    write(out / 'RESULT.json', receipt)
    print(json.dumps({'stage': args.stage, 'status': receipt['status'], 'result': pin(out / 'RESULT.json')}), flush=True)
    return 0 if receipt['status'] in ('complete', 'passed') else 1


if __name__ == '__main__':
    raise SystemExit(main())
