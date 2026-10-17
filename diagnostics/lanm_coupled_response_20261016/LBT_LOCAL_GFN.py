#!/usr/bin/env python3
"""Two missing 7CCN matched-core vacuum GFN2 calls and four-cell comparison."""
from concurrent.futures import ThreadPoolExecutor
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import sys
import time
import traceback

HA = 627.509474
EV = 23.06054783061903
METALS = ('La', 'Dy')
SOURCES = ('7CCO', '7CCN')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(p):
    path = Path(p['path']).resolve(strict=True)
    require(path.is_file() and pin(path)['sha256'] == p['sha256'], 'Changed pin: ' + p['path'])
    return path


def read(p):
    return json.loads(verify(p).read_text())


def write(path, obj):
    with Path(path).open('x') as handle:
        json.dump(obj, handle, indent=2, allow_nan=False)
        handle.write('\n')


def module(name, p):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def command(executable):
    # Exact accepted 12488 default single-point command: no optimization/gradient.
    return [executable, 'core.xyz', '--gfn', '2', '--chrg', '-1', '--uhf', '0',
            '--etemp', '300', '--acc', '0.0001', '--iterations', '500',
            '--norestart', '--json', '--pop']


def environment(task):
    env = os.environ.copy()
    for key, value in task['environment'].items():
        if value is None:
            env.pop(key, None)
        else:
            env[key] = value
    for key in ('OMP_PLACES', 'OMP_PROC_BIND'):
        env.pop(key, None)
    return env


def task_check(task, endpoint, expected):
    require(task['backend'] == 'standalone' and task['medium'] == 'vacuum'
            and task['metal'] == endpoint['metal'] and task['atom_count'] == 48
            and task['charge'] == -1 and task['multiplicity'] == 1 and task['uhf'] == 0
            and task['physical_multiplicity'] == endpoint['multiplicity'], 'GFN physical/effective state differs')
    require(task['source_xyz']['sha256'] == endpoint['xyz']['sha256'] == task['xyz']['sha256'],
            'GFN source coordinates differ')
    for p in (task['source_xyz'], task['xyz'], task['executable'], task['parameters'], *task['runtime_libraries']):
        verify(p)
    require(task['command'] == command(task['executable']['path']), 'GFN command differs')
    require(expected['atoms'] == 48 and expected['charge'] == -1 and expected['multiplicity'] == 1
            and expected['valence_electrons'] % 2 == 0, 'GFN valence state differs')
    require(all(task['environment'].get(k) == '1' for k in
            ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')),
            'Each endpoint requires one thread')


def admit(directory, task, execution, expected, backend, tolerances):
    require(execution['returncode'] == 0 and execution['owned_process_group_empty']
            and execution['command'] == task['command'], 'Native process failed')
    verify(task['xyz'])
    observation = backend.parse_backend(directory, 'xtb', task)
    require(observation['status'] == 'available', observation.get('failure_reason', 'Parser unavailable'))
    payload = json.loads((directory / 'xtbout.json').read_text())
    energy, charges = observation['energy_hartree'], observation['charges_e']
    require(payload['number of electrons'] == expected['valence_electrons']
            and payload['number of unpaired electrons'] == 0
            and payload['method'] == 'GFN2-xTB' and payload['xtb version'].startswith('6.7.1'), 'Native model/state differs')
    require(abs(payload['total energy'] - energy) <= tolerances['xtb_json_stdout_hartree'], 'Native energy disagreement')
    text = (directory / 'stdout.txt').read_text()
    combined = text + '\n' + (directory / 'stderr.txt').read_text()
    require('normal termination of xtb' in combined.lower()
            and not re.search(r'failed to converge|NOT CONVERGED|abnormal termination', combined, re.I), 'Native convergence failed')
    for name, wanted in ((r'net charge', -1), (r'unpaired electrons', 0),
                         (r'# electrons', expected['valence_electrons']),
                         (r'electronic temp\.', 300), (r':  accuracy', .0001)):
        values = re.findall(name + r'\s+([-+\d.Ee]+)', text)
        require(bool(values) and float(values[-1]) == wanted, 'Native setting differs: ' + name)
    for pattern in (r'max\. iterations\s+500', r'restarted\?\s+false',
                    r'PC potential\s+false', r'convergence criteria satisfied after'):
        require(re.search(pattern, text), 'Native setting differs: ' + pattern)
    require('Solvation model:               ALPB' not in text, 'Unexpected solvent')
    require(math.isfinite(energy) and len(charges) == 48 and all(math.isfinite(q) for q in charges), 'Invalid energy/populations')
    closure = abs(math.fsum(charges) + 1)
    require(closure <= tolerances['charge_closure_e'], 'Charge closure failed')
    return {'status': 'accepted', 'energy_hartree': energy, 'charges_e': charges,
            'charge_closure_e': closure, 'standalone_observation': observation,
            'standalone_json': pin(directory / 'xtbout.json'), 'charges': pin(directory / 'charges')}


def inputs(cfg):
    verify(cfg['implementation']); verify(cfg['plan']); verify(cfg['python'])
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Implementation differs')
    require(pin(sys.executable)['sha256'] == cfg['python']['sha256'], 'Python differs')
    compare = module('local_reference_preparation', cfg['preparation_helper'])
    q = read(cfg['qualification']); qc = read(q['config']); manifest = read(q['manifest'])
    require(q['status'] == 'passed' and q['slurm_job_id'] == '12488' and q['accepted'] == 8,
            'Qualification unavailable')
    for p in (q['implementation'], qc['implementation'], qc['stager'], qc['plan'], *qc['pins'], manifest['profile']):
        verify(p)
    require(q['implementation']['sha256'] == qc['implementation']['sha256'], 'Qualification implementation differs')
    require(qc['tolerances']['charge_closure_e'] == 1e-5
            and qc['tolerances']['xtb_json_stdout_hartree'] == 1e-7, 'Original tolerances differ')
    backend = module('local_reference_backend', qc['backend_helper'])
    mace = read(cfg['MACE_result']); mc = read(mace['config']); preparation = read(mc['preparation'])
    require(mace['status'] == 'complete' and mace['slurm_job_id'] == '12471'
            and mace['model_parameters_unchanged'] and mace['model'] == mc['model']
            and mace['versions'] == mc['expected_versions'], 'MACE receipt/model differs')
    for p in (mc['implementation'], mc['plan'], mc['python'], mc['model']['checkpoint'], *mc['helper_pins'].values()):
        verify(p)
    model = mc['model']
    require(model['family'] == 'MACE_OMOL' and model['head'] == 'omol' and model['dtype'] == 'float64'
            and model['solvent'] is None and not model['extra_dispersion'] and not model['periodic']
            and model['execution_adapter'] is None and model['spin_input_convention'] == 'multiplicity', 'MACE method differs')
    pools = compare.preparations(preparation, mc)
    require(pools['matched_core']['compatible'], 'Sources have different local chemical inventory')
    profile = read(manifest['profile']); endpoints = {}; mace_rows = {}; old = {}; expected = {}
    for source in SOURCES:
        for metal in METALS:
            matches = [t for t in mc['tasks'] if (t['source_id'], t['representation'], t['metal']) == (source, 'matched_core', metal)]
            require(len(matches) == 1, 'Missing/duplicate prepared local state')
            t = matches[0]; coordinates = compare.xyz(t['xyz'])
            require(len(coordinates) == 48 and t['charge'] == -1 and t['metal_index'] == 0
                    and t['multiplicity'] == (1 if metal == 'La' else 6) and t['oxidation_state'] == 3, 'Local chemistry differs')
            row = [r for r in mace['rows'] if r['task_id'] == t['task_id']]
            require(len(row) == 1 and row[0]['status'] == 'complete' and row[0]['call_attempted']
                    and all(row[0].get(k) == v for k, v in t.items()) and math.isfinite(row[0]['energy_eV']), 'Saved MACE cell differs')
            compare.state_check(row[0], coordinates, -1, t['multiplicity'], metal, 0, t['all_electron_count'])
            endpoints[(source, metal)] = t; mace_rows[(source, metal)] = row[0]
            valence = sum(sum(profile['element'][a[0]]['refocc']) for a in coordinates) + 1
            basis = sum(sum({'s': 1, 'p': 3, 'd': 5, 'f': 7}[sh[-1]] for sh in profile['element'][a[0]]['shells']) for a in coordinates)
            require(int(valence) == valence and valence % 2 == 0, 'Invalid effective electron count')
            expected[(source, metal)] = {'atoms': 48, 'charge': -1, 'multiplicity': 1,
                                         'valence_electrons': int(valence), 'basis_functions': basis}
        require(compare.paired(*(compare.xyz(endpoints[(source, m)]['xyz']) for m in METALS), 0), 'Paired coordinates differ')
    for metal in METALS:
        rows = [r for r in q['rows'] if (r['task']['metal'], r['task']['medium'], r['task']['backend']) == (metal, 'vacuum', 'standalone')]
        require(len(rows) == 1 and rows[0]['status'] == 'accepted', 'Saved 7CCO GFN missing')
        row = rows[0]; t = row['task']; d = Path(t['directory'])
        require(read(pin(d / 'RESULT.json')) == row and row['expected'] == expected[('7CCO', metal)], 'Archived GFN receipt differs')
        require(read(t['task_record']) == {k: v for k, v in t.items() if k != 'task_record'}, 'Archived TASK differs')
        task_check(t, endpoints[('7CCO', metal)], row['expected'])
        parsed = admit(d, t, read(row['execution']), row['expected'], backend, qc['tolerances'])
        require(parsed['energy_hartree'] == row['energy_hartree'] and parsed['charges_e'] == row['charges_e']
                and parsed['standalone_observation'] == row['standalone_observation'], 'Saved GFN replay differs')
        old[metal] = row
    template = old['La']['task']
    keys = ('command', 'environment', 'executable', 'parameters', 'runtime_libraries', 'tolerances')
    require(all(template[k] == old['Dy']['task'][k] for k in keys), 'Metal runtimes/settings differ')
    actual = backend.dynamic_pins(template['executable']['path'], environment(template))
    require(sorted((p['path'], p['sha256']) for p in actual['libraries']) ==
            sorted((str(verify(p)), p['sha256']) for p in template['runtime_libraries']), 'Loaded runtime libraries differ')
    return dict(backend=backend, endpoints=endpoints, MACE=mace_rows, old=old, expected=expected,
                tolerances=qc['tolerances'], runtime=actual, preparation=mc['preparation'], qualification=cfg['qualification'])


def prepare(metal, world, output):
    endpoint = world['endpoints'][('7CCN', metal)]; d = output / ('7CCN_' + metal); d.mkdir()
    (d / 'core.xyz').write_bytes(verify(endpoint['xyz']).read_bytes())
    t = copy.deepcopy(world['old'][metal]['task']); t.pop('task_record')
    t.update(id='7CCN_' + metal + '_vacuum', directory=str(d), xyz=pin(d / 'core.xyz'),
             source_xyz=endpoint['xyz'], source_id='7CCN', physical_source=endpoint['physical_source'],
             source_cap_map=endpoint['source_cap_map'], prepared_endpoint=endpoint,
             target_index=0, runtime_qualified=True, runtime_qualified_by=world['qualification'])
    task_check(t, endpoint, world['expected'][('7CCN', metal)])
    write(d / 'TASK.json', t); t['task_record'] = pin(d / 'TASK.json')
    return {'task': t, 'expected': world['expected'][('7CCN', metal)], 'status': 'unavailable',
            'attempted': False, 'energy_hartree': None}


def execute(row, world):
    t = row['task']; d = Path(t['directory'])
    try:
        row['attempted'] = True
        execution = world['backend'].call(t['command'], d, environment(t))
        row['execution'] = pin(d / 'EXECUTION.json')
        row.update(admit(d, t, execution, row['expected'], world['backend'], world['tolerances']))
        for p in (t['executable'], t['parameters'], *t['runtime_libraries']): verify(p)
    except Exception as exc:
        row.update(status='unavailable', failure=repr(exc), traceback=traceback.format_exc())
    row['artifacts'] = {p.name: pin(p) for p in sorted(d.iterdir()) if p.is_file()}
    write(d / 'RESULT.json', row)
    return row


def run(config_path, output):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm allocation required')
    require(int(os.environ['SLURM_CPUS_PER_TASK']) == 2, 'Exactly two shared CPUs required')
    cfg = json.loads(config_path.read_text()); output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic(); cpu = time.process_time()
    result = {'status': 'unavailable', 'protocol_id': 'lbt3_local_source_response_GFN_v1',
              'implementation': pin(__file__), 'config': pin(config_path), 'plan': cfg['plan'],
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'declared_new_GFN_calls': 2,
              'new_GFN_calls_attempted': 0, 'reused_GFN_cells': 0, 'reused_MACE_cells': 0,
              'new_MACE_calls': 0, 'new_PBE0_calls': 0, 'solvent_calls': 0, 'optimizer_steps': 0,
              'qualification': cfg['qualification'], 'MACE_result': cfg['MACE_result'],
              'rows': [], 'matrix': [], 'source_strains': {}}
    write(output / 'START.json', result)
    try:
        world = inputs(cfg)
        result.update(preflight='passed', runtime=world['runtime'], preparation=world['preparation'],
                      reused_GFN_cells=2, reused_MACE_cells=4)
        # Complete both coordinate/state/receipt preparations before either call.
        staged = [prepare(m, world, output) for m in METALS]
        with ThreadPoolExecutor(max_workers=2) as pool:
            result['rows'] = list(pool.map(lambda row: execute(row, world), staged))
        result['new_GFN_calls_attempted'] = sum(r['attempted'] for r in result['rows'])
        new = {r['task']['metal']: r for r in result['rows']}
        for source in SOURCES:
            for metal in METALS:
                row = world['old'][metal] if source == '7CCO' else new[metal]
                result['matrix'].append({'source_id': source, 'metal': metal,
                    'endpoint': world['endpoints'][(source, metal)], 'GFN_status': row['status'],
                    'GFN_energy_hartree': row['energy_hartree'] if row['status'] == 'accepted' else None,
                    'GFN_result': pin(Path(row['task']['directory']) / 'RESULT.json'),
                    'GFN_reused': source == '7CCO', 'MACE_energy_eV': world['MACE'][(source, metal)]['energy_eV'],
                    'MACE_native_state': world['MACE'][(source, metal)]['input_state_check']})
        for method, key, factor in (('MACE', 'MACE_energy_eV', EV), ('GFN2', 'GFN_energy_hartree', HA)):
            energies = {(r['source_id'], r['metal']): r[key] for r in result['matrix']}
            strains = {m: (math.fsum([energies[('7CCN', m)], -energies[('7CCO', m)]]) * factor
                           if all(energies[(s, m)] is not None for s in SOURCES) else None) for m in METALS}
            differential = (math.fsum([energies[('7CCN', 'Dy')], -energies[('7CCO', 'Dy')],
                                      -energies[('7CCN', 'La')], energies[('7CCO', 'La')]]) * factor
                            if all(v is not None for v in energies.values()) else None)
            result['source_strains'][method] = {'definition': 'E(7CCN)-E(7CCO), same metal',
                'La_kcal_mol': strains['La'], 'Dy_kcal_mol': strains['Dy'], 'Dy_minus_La_kcal_mol': differential}
        result['status'] = 'complete' if all(r['status'] == 'accepted' for r in result['rows']) else 'unavailable'
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-start, process_CPU_seconds=time.process_time()-cpu,
                  interpretation='Local fixed geometry electronic source response. PBE0 pending separately; no affinity or physical Dy-state validation.',
                  interpretation_resolution_kcal_mol=1., effective_GFN_multiplicity=1, physical_multiplicities={'La': 1, 'Dy': 6})
    write(output / 'RESULT.json', result)
    lines = ['# Local four-cell MACE/GFN source comparison', '',
             'S = E(7CCN) − E(7CCO), same-metal native differences; kcal/mol.', '',
             '| Method | La S | Dy S | Dy S − La S |', '|---|---:|---:|---:|']
    for method, s in result['source_strains'].items():
        lines.append(f"| {method} | {s['La_kcal_mol']} | {s['Dy_kcal_mol']} | {s['Dy_minus_La_kcal_mol']} |")
    lines += ['', f"Status: {result['status']}. New GFN calls attempted: {result['new_GFN_calls_attempted']}/2.",
              f"Reused MACE cells: {result['reused_MACE_cells']}; reused 7CCO GFN cells: {result['reused_GFN_cells']}.",
              'No solvent, optimization or new MACE/PBE0.',
              'PBE0 comparison remains separate. Effective GFN singlet differs from physical Dy sextet.']
    (output / 'REPORT.md').write_text('\n'.join(lines) + '\n')
    return result['status'] == 'complete'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(0 if run(args.config.resolve(), args.output.resolve()) else 1)
