"""Finite whole-peptide/matched-core La/Dy native OMOL scout; Slurm only.

Uses the same native loading/call/state path as lanm_series_followup.mace.
No optimizer, solvent, charge adapter or model modification.
"""
import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import resource
import socket
import sys
import time
import traceback


def pin(path):
    p = Path(path).resolve()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return {'path': str(p), 'sha256': h.hexdigest()}


def verify(record):
    if pin(record['path']) != record:
        raise ValueError('Artifact changed: ' + record['path'])
    return Path(record['path'])


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def run(config_path, out):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'), 'Allocated GPU required')
    c = json.loads(config_path.read_text())
    require(pin(__file__) == c['implementation'], 'Wrong runner')
    verify(c['plan']); verify(c['preparation']); verify(c['native_precedent'])
    for p in c['helper_pins'].values():
        verify(p)
    require(c['model']['execution_adapter'] is None and c['model']['solvent'] is None,
            'Native vacuum model only')
    require(c['model']['dtype'] == 'float64' and not c['model']['extra_dispersion'], 'Model recipe changed')
    tasks = c['tasks']
    expected = {(s, r, m) for s in ('7CCO', '7CCN') for r in ('whole', 'matched_core') for m in ('La', 'Dy')}
    require(len(tasks) == 8 and {(t['source_id'], t['representation'], t['metal']) for t in tasks} == expected,
            'Eight-cell denominator changed')
    require(len({t['task_id'] for t in tasks}) == len(tasks), 'Duplicate task IDs')
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    receipt = {'status': 'running', 'config': pin(config_path), 'plan': c['plan'],
               'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
               'started_epoch': started, 'new_model_calls_attempted': 0,
               'new_solvent_calls': 0, 'new_QM_calls': 0, 'optimizer_steps': 0,
               'declared_cells': 8, 'classifier': None, 'affinity': None}
    write(out / 'START.json', receipt)
    rows = []
    try:
        import torch
        from ase import Atoms
        from ase.data import atomic_numbers
        from mace.calculators import mace_omol
        sys.path.insert(0, str(verify(c['helper_pins']['mace_omol.py']).parent))
        from mace_omol import input_batch
        versions = {k: importlib.metadata.version(k) for k in ('torch', 'mace-torch', 'ase')}
        require(versions == c['expected_versions'], 'Runtime differs from native precedent')
        require(pin(sys.executable) == c['python'], 'Unexpected Python executable')
        require(torch.cuda.is_available() and torch.cuda.device_count() == 1, 'Expected one usable allocated GPU')
        torch.set_num_threads(int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
        torch.set_default_dtype(torch.float64)
        torch.cuda.reset_peak_memory_stats()
        coordinates = {}
        for t in tasks:
            text = verify(t['xyz']).read_text().splitlines()
            atoms = [line.split() for line in text[2:] if line.strip()]
            require(int(text[0]) == len(atoms) == t['atom_count'] and all(len(a) == 4 for a in atoms), 'XYZ size invalid')
            syms = [a[0] for a in atoms]
            q = [[float(x) for x in a[1:]] for a in atoms]
            require(all(math.isfinite(x) for row in q for x in row), 'Nonfinite XYZ')
            mi = t['metal_index']
            require(syms[mi] == t['metal'] and syms.count(t['metal']) == 1, 'Target identity/count invalid')
            require(set(syms) <= {'C', 'H', 'N', 'O', t['metal']}, 'Unexpected molecule element')
            require(t['charge'] == (-2 if t['representation'] == 'whole' else -1), 'Representation charge differs')
            require(t['multiplicity'] == (1 if t['metal'] == 'La' else 6), 'Physical spin hypothesis differs')
            ne = sum(atomic_numbers[s] for s in syms) - t['charge']
            require((ne - t['multiplicity'] + 1) % 2 == 0, 'All-electron spin parity invalid')
            coordinates[(t['source_id'], t['representation'], t['metal'])] = (syms, q, mi)
        for source in ('7CCO', '7CCN'):
            for representation in ('whole', 'matched_core'):
                a, b = [coordinates[(source, representation, m)] for m in ('La', 'Dy')]
                require(a[1] == b[1] and a[2] == b[2] and
                        a[0][:a[2]] + a[0][a[2]+1:] == b[0][:b[2]] + b[0][b[2]+1:], 'Paired coordinate/inventory mismatch')
        tick = time.monotonic()
        calc = mace_omol(model=str(verify(c['model']['checkpoint'])), device='cuda', default_dtype='float64')
        model = calc.models[0]
        require(type(model).__name__ == 'ScaleShiftMACE' and dict(model.embedding_specs) == c['model']['embedding_specs']
                and list(model.heads) == ['omol'] and calc.head == 'omol' and calc.energy_units_to_eV == 1.,
                'Native checkpoint/head/embeddings/units differ')
        for p in model.parameters():
            p.requires_grad_(False)
        receipt['model_load_seconds'] = time.monotonic() - tick
        versions_before = [p._version for p in model.parameters()]
        stop_reason = None
        for t in tasks:
            row = dict(t, status='unavailable', energy_eV=None, call_attempted=False)
            tick = time.monotonic()
            if stop_reason:
                row['reason'] = 'not attempted after preceding failure: ' + stop_reason
            else:
                try:
                    syms, q, mi = coordinates[(t['source_id'], t['representation'], t['metal'])]
                    require({atomic_numbers[s] for s in syms} <= set(model.atomic_numbers.tolist()), 'Element absent from model')
                    atoms = Atoms(syms, positions=q, pbc=False)
                    atoms.info.update(charge=t['charge'], spin=t['multiplicity'])
                    batch = calc._atoms_to_batch(atoms)
                    row['input_state_check'] = input_batch(calc, atoms, t['charge'], t['multiplicity'], batch, mi)
                    with torch.no_grad():
                        row['call_attempted'] = True
                        receipt['new_model_calls_attempted'] += 1
                        v = model(batch.to_dict(), training=False, compute_force=False,
                                  compute_virials=False, compute_stress=False, compute_displacement=False,
                                  compute_hessian=False, compute_edge_forces=False, compute_atomic_stresses=False)
                    energy = float(v['energy'].detach().cpu().item())
                    require(math.isfinite(energy) and v['forces'] is None, 'Nonfinite energy/unexpected forces')
                    require([p._version for p in model.parameters()] == versions_before, 'Model parameters changed')
                    row.update(status='complete', energy_eV=energy)
                except Exception as e:
                    stop_reason = f'{type(e).__name__}: {e}'
                    row.update(status='failed', reason=stop_reason, traceback=traceback.format_exc())
            row['wall_seconds'] = time.monotonic() - tick
            rows.append(row)
            write(out / (t['task_id'] + '.json'), row)
        receipt.update(status='complete' if all(r['status'] == 'complete' for r in rows) else 'partial',
                       rows=rows, model=c['model'], peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated(),
                       device=torch.cuda.get_device_name(),
                       versions=versions,
                       torch_threads=torch.get_num_threads(), model_parameters_unchanged=True)
        require(pin(config_path) == receipt['config'], 'Config changed during execution')
    except Exception as e:
        present = {r['task_id'] for r in rows}
        rows.extend(dict(t, status='unavailable', energy_eV=None, call_attempted=False,
                         reason='run setup/collection failed: ' + str(e))
                    for t in tasks if t['task_id'] not in present)
        receipt.update(status='failed', rows=rows, reason=f'{type(e).__name__}: {e}', traceback=traceback.format_exc())
    receipt.update(finished_epoch=time.time(), peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   elapsed_seconds=time.time() - started,
                   allocation={k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK',
                               'SLURM_MEM_PER_CPU', 'SLURM_MEM_PER_NODE', 'SLURM_GPUS_ON_NODE', 'CUDA_VISIBLE_DEVICES')})
    write(out / 'RESULT.json', receipt)
    return 0 if receipt['status'] == 'complete' else 1


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    raise SystemExit(run(a.config.resolve(), a.output.resolve()))
