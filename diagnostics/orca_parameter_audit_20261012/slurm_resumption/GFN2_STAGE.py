#!/usr/bin/env python3
"""Slurm placement only; use the frozen Ca-transfer solvent implementation."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import queue
import socket
import time

ROOT = Path(__file__).resolve().parents[3]
DRIVER = ROOT / 'diagnostics/orca_parameter_audit_20261012/ca_conditioned_execution/GFN2.py'
DRIVER_SHA256 = 'c8c92032ca0e08541f8e4501da0d92a724b492f297d773d2c6dfbc3083d26a20'
DEPENDENCIES = {
    'scripts/context_factorial_gfn2.py': '4c534bf7ac264c8babd1e1c8aa28c88b8932bfc53708a52666772339b9dcb95c',
    'diagnostics/orca_parameter_audit_20261012/expanded_transfer_inventory/REEVALUATE.py':
        '5f4424023dedb66b50522b5b9b7689b8883e9429356012f837a61a4f13cea152',
}
THREAD_ENV = {
    'OMP_NUM_THREADS': '1', 'OMP_THREAD_LIMIT': '1', 'OMP_DYNAMIC': 'FALSE',
    'OMP_MAX_ACTIVE_LEVELS': '1', 'OPENBLAS_NUM_THREADS': '1',
    'MKL_NUM_THREADS': '1', 'MKL_DYNAMIC': 'FALSE', 'BLIS_NUM_THREADS': '1',
    'NUMEXPR_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1',
}


def allocation():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('This stage requires a real Slurm allocation')
    cpus = int(os.environ['SLURM_CPUS_ON_NODE'])
    nodes = int(os.environ.get('SLURM_JOB_NUM_NODES', os.environ.get('SLURM_NNODES', '1')))
    affinity = sorted(os.sched_getaffinity(0))
    if nodes != 1 or cpus != 12 or len(affinity) < cpus:
        raise RuntimeError(f'Unsupported allocation: nodes={nodes}, CPUs={cpus}, affinity={affinity}')
    return {
        'job_id': os.environ['SLURM_JOB_ID'], 'partition': os.environ.get('SLURM_JOB_PARTITION'),
        'node_list': os.environ.get('SLURM_JOB_NODELIST'), 'host': socket.gethostname(),
        'allocated_cpus_on_node': cpus, 'parent_affinity': affinity,
        'worker_cpu_ids': affinity[:cpus], 'threads_per_child': 1,
        'thread_environment': THREAD_ENV,
    }


def driver():
    digest = hashlib.sha256(DRIVER.read_bytes()).hexdigest()
    if digest != DRIVER_SHA256:
        raise RuntimeError('Frozen scientific driver changed; review before executing')
    for relative, expected in DEPENDENCIES.items():
        if hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'Frozen scientific helper changed: {relative}')
    spec = importlib.util.spec_from_file_location('frozen_ca_gfn2_transfer', DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(scope, collection, base):
    placement = allocation()
    os.environ.update(THREAD_ENV)
    g = driver()
    base = base.resolve()
    base.parent.mkdir(parents=True, exist_ok=True)
    # Also cover preparation: a second copy must not create the same manifest.
    with (base.parent / ('.' + base.name + '.SLURM_STAGE.lock')).open('a') as stage_lock:
        fcntl.flock(stage_lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        stage_started = time.time()
        if not (base / 'manifest.json').exists():
            g.prepare(scope, collection, base)
        with (base / 'EXECUTION.lock').open('a') as execution_lock:
            fcntl.flock(execution_lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            g.dry_run(base)
            manifest = g.f.read(base / 'manifest.json')
            if manifest['scope'] != g.f.pin(scope) or manifest['MACE_collection'] != g.f.pin(collection):
                raise RuntimeError('Existing stage belongs to different scope or MACE collection')
            if manifest['execution_disabled']:
                raise RuntimeError('Archived parser fixtures cannot launch science')
            runtime = g.f.read(manifest['runtime'])
            pending, completed, interrupted = g.r.pending_tasks(
                base, [task for task in manifest['tasks'] if task['reuse'] is None])
            attempt = 1 + len(list(base.glob('SLURM_ATTEMPT_*_START.json')))
            stem = f'SLURM_ATTEMPT_{attempt:03d}'
            workers = min(placement['allocated_cpus_on_node'], len(pending))
            receipt = {
                'started_epoch': stage_started, 'placement': placement,
                'concurrent_workers': workers, 'pending_inputs': len(pending),
                'completed_inputs_reused': completed, 'interrupted_unavailable_inputs': interrupted,
                'manifest': g.f.pin(base / 'manifest.json'), 'wrapper': g.f.pin(__file__),
                'scientific_driver': g.f.pin(DRIVER), 'scientific_settings_changed': False,
                'raw_backend_calls_per_input': 2, 'automatic_retries': False,
            }
            g.f.write(base / (stem + '_START.json'), receipt)
            # Preserve a complete denominator even if the subsequent allocation dies.
            g.collect(base, stem + '_PREEXECUTION_COLLECTION.json')
            cpu_queue = queue.SimpleQueue()
            for cpu in placement['worker_cpu_ids'][:workers]:
                cpu_queue.put(cpu)

            def bind_worker():
                # Linux affinity is per-thread; each molecular child inherits this CPU.
                os.sched_setaffinity(0, {cpu_queue.get_nowait()})

            def execute_one(task):
                try:
                    if len(os.sched_getaffinity(0)) != 1:
                        raise RuntimeError('Worker affinity is not exactly one CPU')
                    g.f.execute_point(task, manifest, base, runtime)
                except Exception as exc:
                    return {'id': task['id'], 'reason': str(exc)}
                return None

            errors = []
            interruption = None
            try:
                if workers:
                    with ThreadPoolExecutor(max_workers=workers, initializer=bind_worker) as pool:
                        errors = [error for error in pool.map(execute_one, pending) if error is not None]
            except BaseException as exc:
                interruption = {'type': type(exc).__name__, 'reason': str(exc)}
                raise
            finally:
                finished = time.time()
                g.f.write(base / (stem + '_FINISH.json'), {
                    'finished_epoch': finished, 'elapsed_seconds_before_final_collection': finished - stage_started,
                    'allocated_CPU_seconds_before_final_collection':
                        (finished - stage_started) * placement['allocated_cpus_on_node'],
                    'exceptions': errors, 'interruption': interruption,
                    'note': 'Allocated CPU-seconds are not measured utilization; backend CPU/RSS are in collection receipts.',
                })
                name = 'FINAL_COLLECTION.json' if not (base / 'FINAL_COLLECTION.json').exists() else stem + '_COLLECTION.json'
                result = g.collect(base, name)
                collected = time.time()
                g.f.write(base / (stem + '_COLLECTION_RECEIPT.json'), {
                    'collection': g.f.pin(base / name), 'manifest': g.f.pin(base / 'manifest.json'),
                    'finished_epoch': collected, 'stage_elapsed_seconds': collected - stage_started,
                    'allocated_CPU_seconds_through_collection':
                        (collected - stage_started) * placement['allocated_cpus_on_node'],
                    'scientific_acceptance': 'read per-source outcomes and unavailable cells',
                    'cost': result['cost'], 'placement': placement,
                })
            if errors:
                raise RuntimeError(f'{len(errors)} input exceptions; see {stem}_FINISH.json')
            return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scope', type=Path, required=True)
    parser.add_argument('--collection', type=Path, required=True)
    parser.add_argument('--directory', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.scope.resolve(), args.collection.resolve(), args.directory)
    print(json.dumps({'declared_sources': result['declared_sources'], 'counts': result['counts']}))
