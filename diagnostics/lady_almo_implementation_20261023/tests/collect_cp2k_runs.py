#!/usr/bin/env python3
"""Join terminal CP2K runner receipts to the unchanged saved-only comparator.

Root may schedule this with --dependency=afterany:<coarse_job>:<fine_job>.
No polling, native execution, retry, or inferred process return codes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent
COMPARATOR_SHA256 = 'd48b05cd01506b1fa9cc4fee0974ec57cd582cace833cd2683ba0be3bb152d2a'
RUNNER_SHA256 = 'efc53aeb34c54513c4c525e5f073de648ed92d131b435213ecb29f363744b5e1'


def need(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return {'path': str(path), 'sha256': h.hexdigest()}


def verify(record, checked):
    key = (record['path'], record['sha256'])
    if key not in checked:
        need(pin(record['path']) == record, 'Missing/changed pinned artifact: ' + record['path'])
        checked.add(key)
    return Path(record['path'])


def terminal(path, cutoff, checked):
    need(path.is_file(), 'Terminal runner receipt is absent: ' + str(path))
    receipt_pin = pin(path)
    row = json.loads(path.read_text())
    need(row['status'] in ('native_complete_needs_comparison', 'native_failed', 'failed'),
         'Runner is not terminal: ' + str(path))
    need('elapsed_seconds' in row and isinstance(row.get('artifacts'), list),
         'Runner terminal finalization is incomplete')
    need(row['cutoff_Ry'] == cutoff and row['native_calls'] == 1,
         'Wrong cutoff or no recorded native attempt')
    need(type(row.get('returncode')) is int, 'No actual native return code; collection cannot infer one')
    need(row['qos'] == 'standard_expand' and row['mpi_ranks'] == 1 and row['omp_threads'] == 8,
         'Changed native resource identity')
    need(row['implementation']['sha256'] == RUNNER_SHA256, 'Unrecognized native runner')
    verify(row['implementation'], checked)
    manifest = json.loads(verify(row['manifest'], checked).read_text())
    need(row['manifest']['sha256'] == pin(HERE/'cp2k/MANIFEST.json')['sha256'],
         'Runner uses a different benchmark manifest')
    runtime = json.loads(verify(row['runtime'], checked).read_text())
    need(row['runtime']['sha256'] == pin(HERE/'cp2k/RUNTIME.json')['sha256'],
         'Runner uses a different binary/container contract')
    for name in ('binary', 'image'):
        verify(runtime[name], checked)
    artifacts = {}
    for record in row['artifacts']:
        p = verify(record, checked)
        need(p.parent.resolve() == path.parent.resolve() and p.name not in artifacts,
             'Runner artifact directory or uniqueness differs')
        artifacts[p.name] = record
    input_path = path.parent/'input.inp'
    output_path = path.parent/'endpoint.out'
    need(row['input'] == artifacts.get('input.inp') == pin(input_path), 'Runner input linkage differs')
    for actual, expected in [('input.inp', f'water_pbe_{cutoff}.inp'),
                             ('BASIS_BENCHMARK', 'BASIS_BENCHMARK'),
                             ('POTENTIAL_BENCHMARK', 'POTENTIAL_BENCHMARK')]:
        need(artifacts[actual]['sha256'] == manifest['inputs'][expected]['sha256'],
             'Executed input/asset differs: ' + actual)
    if row['returncode'] == 0:
        need(row['status'] == 'native_complete_needs_comparison' and not row.get('timed_out')
             and row.get('error') is None, 'Zero native code does not have successful runner admission')
        need(artifacts.get('endpoint.out') == pin(output_path), 'Missing/unbound successful native output')
    else:
        need(row['status'] != 'native_complete_needs_comparison', 'Failed native code marked successful')
        # The existing comparator rejects this actual nonzero code before
        # reading energies, even if native startup produced no endpoint.out.
    need(pin(path) == receipt_pin, 'Runner receipt changed during collection')
    return {'receipt': receipt_pin, 'status': row['status'], 'job_id': row['job_id'],
            'returncode': row['returncode'], 'input': str(input_path), 'output': str(output_path),
            'timed_out': bool(row.get('timed_out')), 'native_calls': row['native_calls']}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--coarse-run', type=Path, required=True)
    ap.add_argument('--fine-run', type=Path, required=True)
    ap.add_argument('--pyscf', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    need(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Saved comparison requires actual Slurm')
    need(os.environ.get('SLURM_JOB_QOS') == 'standard_expand', 'Require standard_expand')
    output = a.output.resolve()
    receipt = output.with_name(output.name + '.collection.json')
    need(not output.exists() and not receipt.exists(), 'Refuse to overwrite comparison or receipt')
    output.parent.mkdir(parents=True, exist_ok=True)
    result = {'status': 'failed_collection_inputs', 'implementation': pin(__file__),
              'job_id': os.environ['SLURM_JOB_ID'], 'new_native_calls': 0,
              'comparator_invoked': False, 'runs': {}, 'error': None}
    code = 1
    try:
        comparator = HERE/'compare_cp2k.py'
        result['comparator'] = pin(comparator)
        need(result['comparator']['sha256'] == COMPARATOR_SHA256, 'Comparator source changed')
        checked = set()
        for name, cutoff in [('coarse', 400), ('fine', 600)]:
            result['runs'][name] = terminal(getattr(a, name+'_run').resolve(), cutoff, checked)
        need(result['runs']['coarse']['job_id'] != result['runs']['fine']['job_id'],
             'Two resolutions must have distinct native parent jobs')
        result['pyscf'] = pin(a.pyscf)
        command = [sys.executable, str(comparator), '--pyscf', str(a.pyscf.resolve()),
                   '--output', str(output)]
        for name in ('coarse', 'fine'):
            row = result['runs'][name]
            command += ['--'+name, row['output'], '--'+name+'-input', row['input'],
                        '--'+name+'-exit-code', str(row['returncode'])]
        result['command'] = command
        result['comparator_invoked'] = True
        code = subprocess.run(command, check=False).returncode
        result['comparator_returncode'] = code
        if output.exists():
            result['comparison'] = pin(output)
            result['comparison_status'] = json.loads(output.read_text())['status']
        result['status'] = 'passed_saved_comparison' if code == 0 else 'failed_saved_comparison'
    except Exception:
        result['error'] = traceback.format_exc()
    with receipt.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return 0 if result['status'] == 'passed_saved_comparison' else 1


if __name__ == '__main__':
    sys.exit(main())
