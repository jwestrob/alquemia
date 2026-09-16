"""Task-owned completion collection and explicitly authorized memory-only recovery."""
import argparse
import fcntl
import json
from pathlib import Path
import subprocess
import sys
import time


def new(path, data):
    with Path(path).open('x') as f:
        json.dump(data, f, indent=2, allow_nan=False); f.write('\n')


def classify_failure(collection, accounting):
    logs = []
    for attempt in collection['attempts']:
        if attempt['accepted'] or not attempt.get('receipt'):
            continue
        receipt = json.loads(Path(attempt['receipt']['path']).read_text())
        if receipt['slurm_job_id'] != accounting['job_id']:
            continue
        logs.append(Path(receipt['log']['path']).read_text())
    text = '\n'.join(logs)
    if 'CUDA out of memory' in text or 'torch.cuda.OutOfMemoryError' in text:
        return 'gpu_memory'
    if 'OUT_OF_MEMORY' in accounting['sacct'] or 'std::bad_alloc' in text or 'MemoryError' in text:
        return 'host_memory'
    return 'other_failure_requires_inspection'


def run(root, manifest, first_job):
    root, manifest = Path(root).resolve(), Path(manifest).resolve()
    data = json.loads(manifest.read_text())
    software = json.loads(Path(data['software']['path']).read_text())
    python = software['python']['path']
    worker = data['implementation']['mace_hybrid.py']['path']
    owned = manifest.parent/'continuation'; owned.mkdir(exist_ok=True)
    with (owned/'watch.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        job, shares, mode, generation = first_job, 1, 'native', 1
        while True:
            accounting_path = owned/f'job_{job}_accounting.json'
            subprocess.run([python, str(root/'scripts/affordable_watch.py'), '--job', job,
                            '--output', str(accounting_path)], check=True)
            accounting = json.loads(accounting_path.read_text())
            collection_path = owned/f'job_{job}_collection.json'
            subprocess.run([python, worker, 'collect', '--manifest', str(manifest),
                            '--output', str(collection_path)], check=True)
            collection = json.loads(collection_path.read_text())
            if collection['status'] == 'complete':
                new(owned/'completion.json', {'status': 'all_12_calls_complete', 'last_job': job,
                    'collection': str(collection_path), 'finished_unix': time.time(),
                    'interpretation': 'Physical checks require review; no calibrated score or automatic promotion.'})
                return
            failure = classify_failure(collection, accounting)
            next_shares, next_mode = shares, mode
            if failure == 'gpu_memory' and mode == 'native':
                next_mode = 'host_offload'
            elif failure == 'host_memory' and shares < 8:
                next_shares = shares * 2
            else:
                new(owned/'completion.json', {'status': 'technical_review_required', 'last_job': job,
                    'failure': failure, 'collection': str(collection_path), 'finished_unix': time.time(),
                    'reason': 'GPU memory is not pooled by reserving more GPUs; no model/physics change or unknown-failure retry.'})
                return
            generation += 1
            cpus = 28 * next_shares
            memory_mib = 2063701 * next_shares // 8
            work = manifest.parent.parent
            cmd = ['sbatch', '--parsable', f'--cpus-per-task={cpus}', f'--gres=gpu:{next_shares}',
                   f'--mem={memory_mib}M', '--time=7-00:00:00', f'--export=ALL,MACE_GPU_SHARES={next_shares}',
                   f'--output={work}/pilot_%j.out', f'--error={work}/pilot_%j.err',
                   str(root/'diagnostics/mace_hybrid_20260916/run_pilot.sbatch'), python,
                   str(manifest), next_mode]
            result = subprocess.run(cmd, capture_output=True, text=True)
            receipt = {'command': cmd, 'returncode': result.returncode, 'stdout': result.stdout,
                       'stderr': result.stderr, 'previous_job': job, 'failure': failure,
                       'gpu_shares': next_shares, 'cpus': cpus, 'host_memory_MiB': memory_mib,
                       'memory_mode': next_mode, 'submitted_unix': time.time(),
                       'time_limit_source': 'cluster standard QOS MaxWall=7 days; no project stopping budget'}
            if result.returncode == 0:
                next_job = result.stdout.strip().split(';')[0]
                receipt['job_id'] = next_job
                receipt['scontrol'] = subprocess.check_output(['scontrol', 'show', 'job', next_job], text=True)
            new(owned/f'launch_{generation:02d}.json', receipt)
            if result.returncode != 0:
                raise RuntimeError('memory-recovery submission failed; receipt retained')
            job, shares, mode = next_job, next_shares, next_mode


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True); p.add_argument('--manifest', required=True); p.add_argument('--first-job', required=True)
    args = p.parse_args()
    run(args.root, args.manifest, args.first_job)
