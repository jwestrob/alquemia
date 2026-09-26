"""Allocation-derived memory wrapper for the unchanged ORCA runtime renderer.

For pinned execution copy this file as render_orca_runtime_input.py and the
original renderer as _base_render_orca_runtime_input.py beside it. Pin both
hashes in preparation. Scientific templates are never edited. All Slurm memory
inputs are interpreted as MiB, converted explicitly to ORCA decimal MB.
"""
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re


def _positive(value, name):
    if value is None or not re.fullmatch(r'[1-9][0-9]*', str(value)):
        raise ValueError(f'{name} must be an explicit positive integer, got {value!r}')
    return int(value)


def allocation_policy(environment=None, *, workers=4, ranks=16):
    """Return bounded concurrency and a 25%-reserved per-rank memory ceiling.

    No host /proc memory, zero-memory whole-node guess, or cpus-per-task as MPI
    slots. Slurm must provide actual allocated CPUs and local task slots.
    """
    env = os.environ if environment is None else environment
    cpus = _positive(env.get('SLURM_CPUS_ON_NODE'), 'SLURM_CPUS_ON_NODE')
    local_slots = env.get('SLURM_NTASKS_PER_NODE') or env.get('SLURM_TASKS_PER_NODE')
    if local_slots is not None:
        match = re.fullmatch(r'([1-9][0-9]*)(?:\(x[1-9][0-9]*\))?', str(local_slots))
        if not match:
            raise ValueError('Heterogeneous/unknown local MPI slots unsupported')
        slots = int(match.group(1))
    else:
        if int(env.get('SLURM_JOB_NUM_NODES', env.get('SLURM_NNODES', '1'))) != 1:
            raise ValueError('Multinode allocation requires explicit local MPI slots')
        slots = _positive(env.get('SLURM_NTASKS'), 'SLURM_NTASKS')
    usable = min(cpus, slots)
    workers = min(_positive(workers, 'workers'), usable)
    ranks = min(_positive(ranks, 'ranks'), usable // workers)
    if env.get('SLURM_MEM_PER_NODE') is not None:
        memory = _positive(env['SLURM_MEM_PER_NODE'], 'SLURM_MEM_PER_NODE')
        source = 'SLURM_MEM_PER_NODE'
    else:
        memory = _positive(env.get('SLURM_MEM_PER_CPU'), 'SLURM_MEM_PER_CPU') * cpus
        source = 'SLURM_MEM_PER_CPU * SLURM_CPUS_ON_NODE'
    raw_mb = memory * (1024**2 / 1000000) * .75 / (workers*ranks)
    maxcore = math.floor(raw_mb / 100) * 100
    if maxcore < 100:
        raise ValueError('Allocation cannot supply even 100 MB per ORCA rank')
    return dict(workers=workers, ranks_per_worker=ranks, allocated_cpus=cpus,
                local_mpi_slots=slots, memory_source=source,
                allocated_memory_MiB=memory, reserve_fraction=.25,
                per_rank_decimal_MB_before_rounding=raw_mb, maxcore_MB=maxcore)


def _base():
    here = Path(__file__).resolve()
    path = here.with_name('_base_render_orca_runtime_input.py')
    if not path.exists():
        path = here.with_name('render_orca_runtime_input.py')
    if path.resolve() == here:
        raise ValueError('Pinned wrapper requires _base_render_orca_runtime_input.py')
    spec = importlib.util.spec_from_file_location('_metal_environment_base_renderer', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


def render_runtime_input(template_path, runtime_path, *, nprocs, sidecar_path=None):
    """Preserve renderer signature/schema/template identity; enrich runtime only.

    The executor must use the recorded worker concurrency (default four).
    METAL_ENV_WORKERS may explicitly select another concurrency before launching
    workers. An overcommitted layout fails rather than changing runner behavior.
    """
    workers = _positive(os.environ.get('METAL_ENV_WORKERS', '4'), 'METAL_ENV_WORKERS')
    nprocs = _positive(nprocs, 'nprocs')
    policy = allocation_policy(workers=workers, ranks=nprocs)
    if policy['workers'] != workers or policy['ranks_per_worker'] != nprocs:
        raise ValueError('Requested executor worker/rank layout exceeds actual allocation')
    template_path, runtime_path = Path(template_path), Path(runtime_path)
    template_text = template_path.read_text()
    if re.search(r'^\s*%maxcore\b', template_text, re.I | re.M):
        raise ValueError('Scientific template already declares MaxCore; ambiguous runtime policy')
    base, base_path = _base()
    record = base.render_runtime_input(template_path, runtime_path, nprocs=nprocs,
                                       sidecar_path=sidecar_path)
    rendered = f"%maxcore {policy['maxcore_MB']}\n" + runtime_path.read_text()
    temporary = runtime_path.with_name(runtime_path.name + '.memory.tmp')
    temporary.write_text(rendered)
    temporary.replace(runtime_path)
    record['runtime_input']['sha256'] = base.sha256_file(runtime_path)
    record['memory_policy'] = policy
    record['memory_policy']['base_renderer_sha256'] = hashlib.sha256(base_path.read_bytes()).hexdigest()
    record['memory_policy']['wrapper_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    sidecar = Path(sidecar_path) if sidecar_path is not None else runtime_path.with_suffix('.json')
    temporary = sidecar.with_name(sidecar.name + '.memory.tmp')
    temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    temporary.replace(sidecar)
    return record
