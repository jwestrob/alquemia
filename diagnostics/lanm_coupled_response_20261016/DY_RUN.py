#!/usr/bin/env python3
"""Run the declared eight real LBT3 parameter controls under Slurm only."""
from concurrent.futures import ThreadPoolExecutor
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import time
import traceback


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def set_environment(env):
    for key, value in env.items():
        if value is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = value


def run(config_path, out):
    assert os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required'
    cfg = json.loads(config_path.read_text())
    root = Path(cfg['repository'])
    sys.path.insert(0, str(root / 'scripts'))
    import native_gfn2_parameters as guard
    # The existing call/parser helper has no engine work at import.
    helper = module('dy_existing_backend_helper', guard.check_pin(cfg['backend_helper']))
    raw = module('dy_existing_native_parser', guard.check_pin(cfg['native_parser']))
    stager = module('dy_stage', guard.check_pin(cfg['stager']))
    assert guard.sha256(__file__) == cfg['implementation']['sha256']
    guard.check_pin(cfg['plan'])
    for p in cfg['pins']:
        guard.check_pin(p)
    cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
    workers = min(4, cpus)
    assert workers == cfg['workers'] and cpus >= workers
    if os.environ.get('SLURM_MEM_PER_NODE'):
        memory_mb = int(os.environ['SLURM_MEM_PER_NODE'])
    else:
        memory_mb = int(os.environ['SLURM_MEM_PER_CPU']) * cpus
    assert memory_mb > 0
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    receipt = {'schema': 'lbt3_dy_parameter_execution.v1', 'config': guard.pin(config_path),
               'plan': cfg['plan'], 'implementation': guard.pin(__file__),
               'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
               'cpus': cpus, 'workers': workers, 'memory_mb': memory_mb,
               'declared_endpoints': 8, 'new_MACE_calls': 0, 'optimizer_steps': 0,
               'status': 'running', 'rows': []}
    write(out / 'START.json', receipt)
    try:
        staging = dict(cfg['staging'])
        staging['native_settings'] = dict(staging['native_settings'],
                                         maxcore_mb=int(0.75 * memory_mb / workers))
        staging['plan'] = cfg['plan']
        staging['execution_config'] = guard.pin(config_path)
        write(out / 'STAGING_CONFIG.json', staging)
        stager.stage(out / 'STAGING_CONFIG.json', out / 'endpoints')
        manifest_path = out / 'endpoints/MANIFEST.json'
        manifest = guard.load_json(manifest_path)
        receipt['manifest'] = guard.pin(manifest_path)
        receipt['maxcore_mb_per_worker'] = staging['native_settings']['maxcore_mb']
        tasks = manifest['tasks']
        assert len(tasks) == 8
        native_environment = next(t['environment'] for t in tasks if t['backend'] == 'native')
        assert all(t['environment'] == native_environment for t in tasks if t['backend'] == 'native')
        # Guard compares os.environ; set once before workers, never mutate it in a thread.
        set_environment(native_environment)
        profile = guard.load_json(manifest['profile']['path'])

        def execute(t):
            d = Path(t['directory'])
            row = {'task': t, 'status': 'unavailable', 'attempted': False,
                   'energy_hartree': None, 'charges_e': None}
            try:
                guard.check_pin(t['task_record']); guard.check_pin(t['xyz'])
                symbols = [x.split()[0] for x in Path(t['xyz']['path']).read_text().splitlines()[2:]]
                valence = sum(sum(profile['element'][s]['refocc']) for s in symbols) - t['charge']
                basis = sum(sum({'s': 1, 'p': 3, 'd': 5, 'f': 7}[sh[-1]]
                                for sh in profile['element'][s]['shells']) for s in symbols)
                assert int(valence) == valence and int(valence) % 2 == 0
                row['expected'] = {'atoms': 48, 'charge': -1, 'multiplicity': 1,
                                   'valence_electrons': int(valence), 'basis_functions': basis}
                env = os.environ.copy()
                for key, value in t['environment'].items():
                    if value is None:
                        env.pop(key, None)
                    else:
                        env[key] = value
                for key in ('OMP_PLACES', 'OMP_PROC_BIND'):
                    env.pop(key, None)
                if t['backend'] == 'native':
                    row['preflight'] = guard.preflight(t['guard']['manifest']['path'])
                else:
                    for p in [t['executable'], t['parameters'], *t['runtime_libraries']]:
                        guard.check_pin(p)
                row['attempted'] = True
                execution = helper.call(t['command'], d, env)
                row['execution'] = guard.pin(d / 'EXECUTION.json')
                assert execution['returncode'] == 0 and execution['owned_process_group_empty']
                guard.check_pin(t['xyz'])
                if t['backend'] == 'native':
                    os.link(d / 'stdout.txt', d / 'endpoint.out')
                    row['postflight'] = guard.postflight(t['guard']['manifest']['path'],
                                                        d / 'endpoint.out', execution['returncode'])
                    data = raw.parse_raw(d)
                    row['native_observation'] = data
                    assert data['normal'] and data['converged'] and not data['printed_gate_failures'], 'native SCF gate'
                    assert data['state'] == row['expected'], 'native atom/electron/basis/state mismatch'
                    assert data['stdout_energy_difference_hartree'] <= cfg['tolerances']['stdout_property_hartree']
                    energy, charges = data['energy_hartree'], data['charges_e']
                else:
                    data = helper.parse_backend(d, 'xtb', t)
                    row['standalone_observation'] = data
                    assert data['status'] == 'available', data.get('failure_reason')
                    energy, charges = data['energy_hartree'], data['charges_e']
                    payload = guard.load_json(d / 'xtbout.json')
                    row['standalone_json'] = guard.pin(d / 'xtbout.json')
                    assert payload['number of electrons'] == int(valence)
                    assert payload['number of unpaired electrons'] == 0
                    assert payload['method'] == 'GFN2-xTB' and payload['xtb version'].startswith('6.7.1')
                    assert abs(payload['total energy'] - energy) <= cfg['tolerances']['xtb_json_stdout_hartree']
                    text = (d / 'stdout.txt').read_text()
                    combined = text + '\n' + (d / 'stderr.txt').read_text()
                    assert 'normal termination of xtb' in combined.lower()
                    assert not re.search(r'failed to converge|NOT CONVERGED|abnormal termination', combined, re.I)
                    for name, expected in ((r'net charge', -1), (r'unpaired electrons', 0),
                                           (r'# electrons', int(valence)),
                                           (r'electronic temp\.', 300), (r':  accuracy', 0.0001)):
                        values = re.findall(name + r'\s+([-+\d.Ee]+)', text)
                        assert values and float(values[-1]) == expected, 'actual standalone setting: ' + name
                    for pattern in (r'max\. iterations\s+500', r'restarted\?\s+false',
                                    r'PC potential\s+false', r'convergence criteria satisfied after'):
                        assert re.search(pattern, text), 'actual standalone setting: ' + pattern
                    if t['medium'] == 'ALPB(Water)':
                        for term in ('Solvation model:               ALPB', 'Reference state                gsolv',
                                     'Solvent                        water', 'Grid points                               230',
                                     'Ion screening                  false'):
                            assert term in text, 'actual standalone solvent setting: ' + term
                    else:
                        assert 'Solvation model:               ALPB' not in text
                assert math.isfinite(energy) and len(charges) == 48 and all(math.isfinite(q) for q in charges)
                assert abs(sum(charges) + 1) <= cfg['tolerances']['charge_closure_e']
                row.update(status='accepted', energy_hartree=energy, charges_e=charges,
                           charge_closure_e=abs(sum(charges) + 1))
            except Exception as exc:
                row['failure'] = repr(exc)
                row['traceback'] = traceback.format_exc()
            row['artifacts'] = {p.name: guard.pin(p) for p in d.iterdir()
                                if p.is_file() and p.suffix in ('.json', '.txt', '.inp', '.xyz', '.out')}
            write(d / 'RESULT.json', row)
            print(json.dumps({'id': t['id'], 'status': row['status'], 'failure': row.get('failure')}), flush=True)
            return row

        with ThreadPoolExecutor(max_workers=workers) as executor:
            rows = list(executor.map(execute, tasks))
        receipt['rows'] = rows
        by = {(r['task']['metal'], r['task']['medium'], r['task']['backend']): r for r in rows}
        parity, components = [], []
        for metal in ('La', 'Dy'):
            for medium in ('vacuum', 'ALPB(Water)'):
                a, b = [by[(metal, medium, e)] for e in ('native', 'standalone')]
                item = {'metal': metal, 'medium': medium, 'passed': False, 'energy_delta_hartree': None,
                        'max_charge_delta_e': None}
                if a['status'] == b['status'] == 'accepted':
                    de = a['energy_hartree'] - b['energy_hartree']
                    dq = max(abs(x-y) for x,y in zip(a['charges_e'], b['charges_e']))
                    item.update(energy_delta_hartree=de, max_charge_delta_e=dq,
                                passed=abs(de) <= cfg['tolerances']['energy_parity_hartree']
                                and dq <= cfg['tolerances']['population_parity_e'])
                parity.append(item)
            item = {'metal': metal, 'passed': False, 'difference_hartree': None}
            four = [by[(metal, medium, e)] for e in ('native','standalone') for medium in ('ALPB(Water)','vacuum')]
            if all(r['status'] == 'accepted' for r in four):
                a, b, c, d = [r['energy_hartree'] for r in four]
                delta = math.fsum([a, -b, -c, d])
                item.update(difference_hartree=delta,
                            passed=abs(delta) <= cfg['tolerances']['energy_parity_hartree'])
            components.append(item)
        receipt.update(parity=parity, solvent_component_parity=components,
                       attempted=sum(r['attempted'] for r in rows),
                       accepted=sum(r['status'] == 'accepted' for r in rows),
                       status='passed' if all(p['passed'] for p in parity + components) else 'failed_qualification')
    except BaseException as exc:
        receipt.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    receipt['elapsed_seconds'] = time.monotonic() - start
    receipt['scope'] = 'fixed real-core parameter/runtime qualification only; no registry/production promotion or affinity claim'
    write(out / 'RESULT.json', receipt)
    print(json.dumps({k: receipt[k] for k in ('status','elapsed_seconds')}), flush=True)
    return receipt['status'] == 'passed'


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    if not run(args.config.resolve(), args.output.resolve()):
        raise SystemExit(1)
