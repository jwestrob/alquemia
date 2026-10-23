#!/usr/bin/env python3
"""Native MACE E/F on an explicitly pinned, finite paired chemistry inventory."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import socket
import sys
import time
import traceback


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': digest}


def verify(record):
    if pin(record['path']) != record:
        raise ValueError('Input changed: ' + record['path'])
    return Path(record['path'])


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def xyz(record):
    lines = verify(record).read_text().splitlines()
    rows = [line.split() for line in lines[2:] if line.strip()]
    if int(lines[0]) != len(rows) or any(len(row) != 4 for row in rows):
        raise ValueError('Invalid XYZ')
    symbols = [row[0] for row in rows]
    positions = [[float(x) for x in row[1:]] for row in rows]
    if not all(math.isfinite(x) for row in positions for x in row):
        raise ValueError('Nonfinite XYZ')
    return symbols, positions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Allocated GPU job required')
    config_pin = pin(args.config)
    cfg = json.loads(Path(args.config).read_text())
    verify(cfg['implementation'])
    if cfg['implementation'] != pin(__file__):
        raise ValueError('Wrong implementation')
    for name in ('specification', 'source_inventory', 'native_helper', 'native_config'):
        verify(cfg[name])
    if len(cfg['tasks']) != cfg['maximum_calls'] or len({t['id'] for t in cfg['tasks']}) != len(cfg['tasks']):
        raise ValueError('Finite task inventory differs')
    prepared = []
    for task in cfg['tasks']:
        symbols, positions = xyz(task['xyz'])
        metal, index = task['metal'], task['metal_index']
        if metal not in ('La', 'Dy') or symbols[index] != metal or symbols.count(metal) != 1:
            raise ValueError('Invalid metal inventory')
        if task['physical_multiplicity'] != (1 if metal == 'La' else 6):
            raise ValueError('Physical multiplicity differs')
        if task['charge'] not in (2, 3) or set(symbols) - {'H', 'C', 'N', 'O', metal}:
            raise ValueError('Unexpected chemistry')
        prepared.append((task, symbols, positions))
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    result = {'status': 'running', 'config': config_pin, 'implementation': pin(__file__),
              'job_id': os.environ['SLURM_JOB_ID'], 'hostname': socket.gethostname(),
              'input_admission': 'passed', 'rows': [], 'native_calls': 0,
              'force_units': 'eV/Angstrom', 'energy_units': 'eV',
              'physical_accuracy_qualified': False, 'affinity_claim': False}
    write(out / 'START.json', result)
    try:
        spec = importlib.util.spec_from_file_location('pinned_native', verify(cfg['native_helper']))
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        native_cfg = json.loads(verify(cfg['native_config']).read_text())
        native = helper.NativeMACE(native_cfg)
        result['device'] = native.device
        result['versions'] = native.versions
        result['model'] = native_cfg['model']
        for task, symbols, positions in prepared:
            row = {'id': task['id'], 'task': task, 'status': 'unavailable'}
            tick = time.monotonic()
            try:
                q, mult = task['charge'], task['physical_multiplicity']
                if (sum(native.Z[s] for s in symbols) - q - mult + 1) % 2:
                    raise ValueError('All-electron physical parity differs')
                if not {native.Z[s] for s in symbols} <= set(native.model.atomic_numbers.tolist()):
                    raise ValueError('Unsupported checkpoint element')
                atoms = native.Atoms(symbols, positions=positions, pbc=False)
                atoms.info.update(charge=q, spin=mult)
                batch = native.calc._atoms_to_batch(atoms)
                state = native.input_batch(native.calc, atoms, q, mult, batch, task['metal_index'])
                if state['graph_count'] != 1 or state['metal_neighbor_edge_count'] <= 0:
                    raise ValueError('Unexpected graph/metal environment')
                with native.torch.enable_grad():
                    result['native_calls'] += 1
                    value = native.model(batch.to_dict(), training=False, compute_force=True,
                                         compute_virials=False, compute_stress=False,
                                         compute_displacement=False, compute_hessian=False,
                                         compute_edge_forces=False, compute_atomic_stresses=False)
                energy = float(value['energy'].detach().cpu().item())
                force = value['forces'].detach().cpu()
                if not math.isfinite(energy) or tuple(force.shape) != (len(symbols), 3) or not bool(native.torch.isfinite(force).all()):
                    raise ValueError('Nonfinite energy/force')
                if [p._version for p in native.model.parameters()] != native.parameter_versions:
                    raise ValueError('Model parameters changed')
                row.update(status='complete', energy_eV=energy, forces_eV_A=force.tolist(),
                           input_state=state, force_convention='negative energy gradient')
            except Exception:
                row['error'] = traceback.format_exc()
            row['elapsed_seconds'] = time.monotonic() - tick
            write(out / (task['id'] + '.json'), row)
            result['rows'].append(row)
        result['status'] = 'complete' if all(r['status'] == 'complete' for r in result['rows']) else 'incomplete'
        verify(config_pin)
    except Exception:
        result.update(status='failed', error=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic() - start
    write(out / 'RESULT.json', result)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    sys.exit(main())
