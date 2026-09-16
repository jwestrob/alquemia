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


def run(root, manifest, first_job, partition='gpu_h200', cpus_per_share=28,
        memory_per_share=200000., nodelist=None):
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
            cpus = cpus_per_share * next_shares
            memory_mib = int(memory_per_share * next_shares)
            actual_memories = [int(r['allocated_host_mem_MiB']) for r in collection['rows'].values()
                               if r.get('slurm_job_id') == job and r.get('allocated_host_mem_MiB')]
            # A killed host-OOM attempt may have no worker result; use the
            # scheduler request in accounting as the minimum known allocation.
            if not actual_memories:
                detail = subprocess.check_output(['sacct', '-j', job, '-n', '-P', '--format=JobID,ReqMem'], text=True)
                for line in detail.splitlines():
                    fields = line.split('|')
                    if fields[0] == job:
                        value = fields[1].rstrip('nc')
                        if value.endswith('M'):
                            actual_memories.append(int(value[:-1]))
                        elif value.endswith('G'):
                            actual_memories.append(int(float(value[:-1])*1024))
            minimum_memory = max(actual_memories or [int(memory_per_share)]) + 1 if failure == 'host_memory' else int(memory_per_share)
            work = manifest.parent.parent
            cmd = ['sbatch', '--parsable', f'--partition={partition}', f'--cpus-per-task={cpus}', f'--gres=gpu:{next_shares}',
                   f'--mem={memory_mib}M', '--time=7-00:00:00',
                   f'--export=ALL,MACE_GPU_SHARES={next_shares},MACE_MIN_MEMORY_MIB={minimum_memory}',
                   f'--output={work}/pilot_%j.out', f'--error={work}/pilot_%j.err',
                   *([f'--nodelist={nodelist}'] if nodelist else []),
                   str(root/'diagnostics/mace_hybrid_20260916/run_pilot.sbatch'), python,
                   str(manifest), next_mode]
            result = subprocess.run(cmd, capture_output=True, text=True)
            receipt = {'command': cmd, 'returncode': result.returncode, 'stdout': result.stdout,
                       'stderr': result.stderr, 'previous_job': job, 'failure': failure,
                       'gpu_shares': next_shares, 'cpus': cpus, 'host_memory_MiB': memory_mib,
                       'minimum_host_memory_MiB': minimum_memory,
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
    p.add_argument('--partition', default='gpu_h200'); p.add_argument('--cpus-per-share', type=int, default=28)
    p.add_argument('--memory-per-share', type=float, default=200000.); p.add_argument('--nodelist')
    args = p.parse_args()
    run(args.root, args.manifest, args.first_job, args.partition, args.cpus_per_share,
        args.memory_per_share, args.nodelist)
