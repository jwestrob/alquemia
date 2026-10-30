#!/usr/bin/env python3
"""Submit one explicitly named immutable development run and arm its completion wake."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

SOURCE = Path(__file__).resolve().parent
ROOT = SOURCE.parents[1]
WORK = ROOT / 'workspaces/lady_almo_implementation_20261023'
PYTHON = WORK / 'venv/bin/python'
THREAD = '01a0a63a-ee36-7483-b726-1f7f5b7f75cc'


def stamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--name', required=True)
    p.add_argument('--cpus', type=int, default=8)
    p.add_argument('--ntasks', type=int, default=1,
                   help='MPI task slots within --cpus total allocated CPUs')
    p.add_argument('--memory-gb', type=int, default=32)
    p.add_argument('--hours', type=int, default=2)
    p.add_argument('--afterany', default='', help='Colon-separated terminal job dependencies')
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    if not a.name.replace('_', '').replace('-', '').isalnum():
        p.error('name must be a single alphanumeric path component')
    if not (1 <= a.cpus <= 24 and 1 <= a.memory_gb <= 256 and 1 <= a.hours <= 12):
        p.error('allocation outside this implementation task bounds')
    if not (1 <= a.ntasks <= a.cpus) or a.cpus % a.ntasks:
        p.error('--ntasks must divide --cpus total allocation')
    threads_per_task = a.cpus // a.ntasks
    if a.afterany and not all(j.isdigit() for j in a.afterany.split(':')):
        p.error('--afterany requires colon-separated job IDs')
    command = a.command[1:] if a.command[:1] == ['--'] else a.command
    if not command or not PYTHON.exists():
        p.error('command and installed environment required')
    task = WORK / a.name
    task.mkdir(parents=True, exist_ok=False)
    snapshot = task / 'source'
    snapshot.mkdir()
    pins = {}
    for f in sorted(SOURCE.rglob('*')):
        if f.is_file() and (f.suffix in {'.py', '.json', '.md', '.inp', '.dat'} or
                            f.name in {'BASIS_BENCHMARK', 'POTENTIAL_BENCHMARK'}):
            if '__pycache__' in f.parts:
                continue
            rel = f.relative_to(SOURCE)
            data = f.read_bytes()
            target = snapshot / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
            digest = hashlib.sha256(data).hexdigest()
            assert hashlib.sha256(f.read_bytes()).hexdigest() == digest, f
            pins[str(rel)] = digest
    replacements = {'@PYTHON@': str(PYTHON), '@SOURCE@': str(snapshot),
                    '@TASK@': str(task), '@ROOT@': str(ROOT), '@CPUS@': str(a.cpus)}
    for key, value in replacements.items():
        command = [part.replace(key, value) for part in command]
    environment_pins = {str(f): hashlib.sha256(f.read_bytes()).hexdigest()
                        for f in (WORK/'ENVIRONMENT.json', WORK/'requirements.lock',
                                  WORK/'PIP_REPORT.json') if f.exists()}
    config = {'created_utc': stamp(), 'command': command, 'source_sha256': pins,
              'environment_sha256': environment_pins,
              'cpus': a.cpus, 'ntasks': a.ntasks, 'cpus_per_task': threads_per_task,
              'memory_GiB': a.memory_gb, 'hours': a.hours,
              'qos': 'standard_expand', 'afterany': a.afterany,
              'task': str(task), 'repo': str(ROOT)}
    write_json(task / 'CONFIG.json', config)
    lines = ['#!/bin/bash', '#SBATCH --job-name=lady-' + a.name,
             '#SBATCH --account=banfield', '#SBATCH --partition=standard',
             '#SBATCH --qos=standard_expand', '#SBATCH --nodes=1',
             f'#SBATCH --ntasks={a.ntasks}', f'#SBATCH --cpus-per-task={threads_per_task}',
             # 2026-10-24: no --mem/--time (per-CPU default memory, QOS walltime); --memory-gb/--hours are recorded only
             '#SBATCH --no-requeue', 'set -euo pipefail',
             'export LADY_REPO_ROOT=' + shlex.quote(str(ROOT)),
             f'export OMP_NUM_THREADS={threads_per_task}', f'export MKL_NUM_THREADS={threads_per_task}',
             'export OPENBLAS_NUM_THREADS=1', 'export OMP_MAX_ACTIVE_LEVELS=1',
             'export PYTHONUNBUFFERED=1', 'export PYTHONPATH=' + shlex.quote(str(snapshot)),
             'exec ' + shlex.join(command), '']
    script = task / 'job.sbatch'
    script.write_text('\n'.join(lines))
    argv = ['sbatch', '--parsable', '--output', str(task / 'slurm-%j.out'),
            '--chdir', str(task), str(script)]
    if a.afterany:
        argv.insert(1, '--dependency=afterany:' + a.afterany)
    result = subprocess.run(argv, capture_output=True, text=True, timeout=60)
    submission = {'submitted_utc': stamp(), 'argv': argv, 'stdout': result.stdout,
                  'stderr': result.stderr, 'returncode': result.returncode,
                  'config_sha256': hashlib.sha256((task/'CONFIG.json').read_bytes()).hexdigest()}
    write_json(task / 'SUBMISSION.json', submission)
    result.check_returncode()
    job = result.stdout.strip().split(';')[0]
    assert job.isdigit(), result.stdout
    submission['job_id'] = job
    watcher = [sys.executable, str(ROOT/'diagnostics/slurm_wakeup_20260927/watch_and_queue.py'),
               '--jobs', job, '--thread', THREAD,
               '--codex', '/home/jwestrob/.nvm/versions/node/v22.23.1/bin/codex',
               '--receipt', str(task/'completion_event.json'),
               '--event', 'lady-almo-' + a.name, '--interval', '20']
    with (task/'watch.log').open('w') as output:
        proc = subprocess.Popen(watcher, stdout=output, stderr=subprocess.STDOUT,
                                start_new_session=True, stdin=subprocess.DEVNULL)
    submission['completion_watcher'] = {'pid': proc.pid, 'argv': watcher}
    write_json(task/'SUBMISSION.json', submission)
    print(json.dumps({'job_id': job, 'task': str(task), 'cpus': a.cpus,
                      'memory_GiB': a.memory_gb, 'watcher_pid': proc.pid}))


if __name__ == '__main__':
    main()
