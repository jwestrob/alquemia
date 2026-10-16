#!/usr/bin/env python3
"""Slurm-only saved-output diagnostic; no molecular calls or admission changes.

CLI --config CONFIG --output FRESH. Config pins implementation, plan,
failed_static_result. Run with the failed static job's actual Python environment.
"""
import argparse
import contextlib
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import os
from pathlib import Path
import socket
import sys
import time
import traceback


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return {'path': str(path), 'sha256': h.hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Changed pin: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def difference(np, saved, replayed):
    a, b = np.asarray(saved, dtype=float), np.asarray(replayed, dtype=float)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(), 'Array shape/finite mismatch')
    d = np.abs(a-b)
    indices = np.flatnonzero(a.ravel() != b.ravel())
    worst = sorted(indices.tolist(), key=lambda i: (-float(d.ravel()[i]), i))[:12]
    return {'shape': list(a.shape), 'exact': saved == replayed,
            'differing_values': int(len(indices)), 'max_absolute_difference': float(d.max(initial=0)),
            'largest_differences': [{'index': [int(k) for k in np.unravel_index(i, a.shape)],
                                    'saved': float(a.ravel()[i]), 'replayed': float(b.ravel()[i]),
                                    'replayed_minus_saved': float((b-a).ravel()[i])} for i in worst]}


def diagnose(cfg, result):
    verify(cfg['plan'])
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Implementation differs')
    failed = read(cfg['failed_static_result'])
    require(failed['status'] == 'execution_failure' and failed['slurm_job_id'] == '12673'
            and failed['failure'] == "ValueError('Saved field raw replay differs: mapped_physical_gradient')"
            and all(failed[k] == 0 for k in ('new_MACE_calls_attempted', 'new_ALPB_calls_attempted', 'new_SCC_calls', 'new_MM_calls')),
            'Expected zero-call failed static parent')
    static_cfg = read(failed['config'])
    verify(failed['implementation'])
    field = read(static_cfg['field_result'])
    require(field['status'] == 'complete' and field['slurm_job_id'] == '12570'
            and field['accepted'] == field['attempted'] == 4, 'Field parent unavailable')
    mace = read(static_cfg['mace_force_result']); mcfg = read(mace['config'])
    require(pin(sys.executable)['sha256'] == mcfg['python']['sha256'], 'Pinned Python binary differs')
    local = module('mex_replay_local', verify(field['implementation']))
    world = local.inputs(verify(field['config']))
    np, mapper = world[5], world[8]
    capture = io.StringIO()
    with contextlib.redirect_stdout(capture):
        np.show_config()
    packages = {}
    for name in ('numpy', 'scipy', 'torch', 'mace-torch'):
        try:
            dist = importlib.metadata.distribution(name)
            packages[name] = {'version': dist.version, 'location': str(dist.locate_file(''))}
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    result['runtime'] = {'sys_executable': sys.executable, 'python_pin': pin(sys.executable),
                         'sys_prefix': sys.prefix, 'python_version': sys.version,
                         'packages': packages, 'numpy_file': pin(np.__file__), 'numpy_config': capture.getvalue(),
                         'mapper': pin(mapper.__file__), 'mapper_numpy_same_module': mapper.np is np,
                         'unit_factor_kcal_mol_A_per_Eh_bohr': mapper.FACTOR,
                         'thread_environment': {k: os.environ.get(k) for k in
                             ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS')}}
    result.update(field_result=static_cfg['field_result'], field_config=field['config'],
                  field_implementation=field['implementation'], failed_static_config=failed['config'])
    keys = ('energy_hartree', 'charges_e', 'qm_gradient', 'pointcharge_gradient', 'mapped_physical_gradient')
    rows = []; result['rows'] = rows
    require(len(field['rows']) == 4 and {(r['task']['metal'], r['task']['kind']) for r in field['rows']}
            == {(m,k) for m in ('La','Dy') for k in ('isolated','embedded')}, 'Field denominator differs')
    for saved in field['rows']:
        require(saved['status'] == 'accepted' and read(saved['result']) == {k:v for k,v in saved.items() if k != 'result'},
                'Saved accepted result content differs')
        for p in saved['artifacts'].values():
            verify(p)
        execution = read(saved['execution'])
        replay = local.admit(Path(saved['execution']['path']).parent, saved['task'], execution, world)
        mapping = read(saved['task']['mapping'])
        direct = mapper.scatter(mapping, saved['qm_gradient']['hartree_per_bohr'],
                                saved['pointcharge_gradient']['hartree_per_bohr'])
        require(direct == replay['mapped_physical_gradient'], 'Saved raw gradients and replay scatter disagree')
        a, b = saved['mapped_physical_gradient'], replay['mapped_physical_gradient']
        require(set(a) == set(b) and set(a['components']) == set(b['components']), 'Mapped schema differs')
        fields = {k: difference(np, a[k], b[k]) for k in a if k != 'components'}
        fields.update({'components.'+k: difference(np, a['components'][k], b['components'][k]) for k in a['components']})
        rows.append({'metal': saved['task']['metal'], 'kind': saved['task']['kind'],
                     'saved_result': saved['result'], 'mapping': saved['task']['mapping'],
                     'raw_artifact_pins_verified': True, 'exact_fields': {k: saved[k] == replay[k] for k in keys},
                     'saved_raw_scatter_equals_raw_replay_scatter': True, 'mapped_differences': fields,
                     'raw_QM_gradient': difference(np, saved['qm_gradient']['hartree_per_bohr'], replay['qm_gradient']['hartree_per_bohr']),
                     'raw_PC_gradient': difference(np, saved['pointcharge_gradient']['hartree_per_bohr'], replay['pointcharge_gradient']['hartree_per_bohr'])})
    result.update(status='complete', raw_fields_exact_all=all(all(r['exact_fields'][k] for k in keys[:-1]) for r in rows),
                  mapped_fields_exact_all=all(r['exact_fields']['mapped_physical_gradient'] for r in rows),
                  interpretation='Diagnostic only: no acceptance tolerance or source artifact changed')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    cfg = json.loads(a.config.read_text()); start = time.monotonic()
    result = {'status': 'running', 'implementation': pin(__file__), 'config': pin(a.config),
              'plan': cfg['plan'], 'failed_static_result': cfg['failed_static_result'],
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
              'new_MACE_calls': 0, 'new_GFN2_calls': 0, 'new_ALPB_calls': 0, 'new_MM_calls': 0}
    try:
        diagnose(cfg, result)
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic()-start
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds']*int(os.environ['SLURM_CPUS_PER_TASK'])
    with (out/'RESULT.json').open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False); f.write('\n')
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
