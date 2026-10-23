#!/usr/bin/env python3
"""Root-only finite submission; default prints a plan without submitting."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parents[3]
WATCH = ROOT / 'diagnostics/slurm_wakeup_20260927/watch_and_queue.py'
THREAD = '01a0a63a-ee36-7483-b726-1f7f5b7f75cc'
CODEX = '/home/jwestrob/.nvm/versions/node/v22.23.1/bin/codex'


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(value):
    assert pin(value['path']) == value, f'Changed input: {value}'
    return Path(value['path'])


def write_new(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--release', required=True, type=Path)
    ap.add_argument('--stage', required=True, choices=['scouts', 'remaining'])
    ap.add_argument('--submit', action='store_true')
    args = ap.parse_args()
    r = json.loads(args.release.read_text())
    panel_path = verify(r['panel'])
    source = verify(r['implementation'])
    review = verify(r['independent_review'])
    p = json.loads(panel_path.read_text())
    wrapper = verify(p['run_wrapper'])
    cells = r['released_cell_ids']
    assert len(p['cells']) == 68 and r['maximum_native_calls'] == 68
    assert p['ranks'] == r['ranks'] == 4
    assert len(cells) == len(set(cells))
    assert set(cells) <= {c['cell_id'] for c in p['cells']}
    scouts = {f'pbe0__aquo9__Dy__water_01__{role}' for role in ['metal_water', 'water_only']}
    if args.stage == 'scouts':
        assert set(cells) == scouts
    else:
        assert len(cells) == 66 and not set(cells) & scouts
    run = Path(r['execution_root'])
    assert run.is_absolute() and run.name == 'run_v1'
    records = []
    for cell in cells:
        task = run / 'tasks' / cell
        assert not (task / 'SUBMISSION_INTENT.json').exists(), f'Prior intent: inspect, never duplicate {task}'
        assert not (task / 'execution').exists(), f'Prior execution: {task}'
        argv = ['sbatch', '--parsable', f'--job-name=cb-{cell}',
                f'--output={task}/slurm-%j.out', f'--error={task}/slurm-%j.err',
                str(wrapper), str(source), 'run', '--panel', str(panel_path),
                '--release', str(args.release.resolve()), '--cell', cell,
                '--out', str(task / 'execution')]
        records.append({'cell_id': cell, 'argv': argv, 'task': str(task),
                        'release': pin(args.release), 'panel': pin(panel_path)})
    print(json.dumps({'stage': args.stage, 'new_jobs': len(records), 'ranks_each': 4,
                      'native_calls_total_bound': 68, 'review': str(review),
                      'first_command': records[0]['argv'], 'execute': args.submit}, indent=2), flush=True)
    if not args.submit:
        return
    assert WATCH.is_file() and Path(CODEX).is_file()
    queue = subprocess.run(['squeue', '-h', '-u', 'jwestrob', '-o', '%i'],
                           text=True, capture_output=True, check=True).stdout
    assert len(queue.splitlines()) + len(records) <= 200, 'Submission cap would be exceeded'
    stage_dir = run / args.stage
    stage_dir.mkdir(parents=True, exist_ok=False)
    submitted = []
    error = None
    try:
        for record in records:
            task = Path(record['task'])
            task.mkdir(parents=True, exist_ok=True)
            record['utc'] = now()
            write_new(task / 'SUBMISSION_INTENT.json', record)
            result = subprocess.run(record['argv'], text=True, capture_output=True)
            record.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
            write_new(task / 'SBATCH_RETURN.json', record)
            if result.returncode != 0:
                raise RuntimeError(f'Submission failed for {record["cell_id"]}: {result.stderr}')
            job = result.stdout.strip().split(';')[0]
            assert job.isdigit(), 'Ambiguous submission result; inspect scheduler before retry'
            record['job_id'] = job
            write_new(task / 'SUBMISSION.json', record)
            submitted.append(record)
            print(f'{record["cell_id"]}: {job}', flush=True)
    except BaseException as exc:
        error = repr(exc)
    finally:
        write_new(stage_dir / 'SUBMISSION.json', {
            'utc': now(), 'stage': args.stage, 'tasks': submitted,
            'jobs': [x['job_id'] for x in submitted], 'error': error,
            'release': pin(args.release), 'submitter': pin(__file__),
            'memory_time': 'scheduler defaults', 'qos': 'standard_expand'})
        if submitted:
            argv = ['/usr/bin/python3', str(WATCH), '--jobs',
                    *[x['job_id'] for x in submitted], '--thread', THREAD,
                    '--codex', CODEX, '--receipt', str(stage_dir / 'completion_event.json'),
                    '--event', f'lady-common-basis-{args.stage}-v1', '--interval', '15']
            with (stage_dir / 'watcher.log').open('x') as log:
                proc = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT,
                                        start_new_session=True)
            write_new(stage_dir / 'WATCHER.json', {'pid': proc.pid, 'argv': argv, 'utc': now()})
    if error:
        raise RuntimeError(error)


if __name__ == '__main__':
    main()
