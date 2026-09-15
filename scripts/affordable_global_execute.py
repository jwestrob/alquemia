"""Memory-aware launch of the approved two-endpoint native whole-chain pilot.

Native startup was observed at 32,055,664 KiB RSS per rank. Reserve 64 GiB
per rank for launch planning, leave 25% node RAM unassigned, and retain the
existing runner's 16-rank endpoint default. These are hardware-fit controls,
not elapsed-time, accumulated-cost, or scientific-iteration budgets.
"""
import argparse
import os
from pathlib import Path
import socket

from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from run_orca_task_manifest import run_manifest, _parse_allocation


def layout(cpus, memory_bytes):
    ranks = min(16, cpus//2, int(memory_bytes*0.75)//(2*64*1024**3))
    if ranks < 1:
        raise InvalidArtifact('node cannot fit two endpoints with measured native per-rank memory allowance')
    return {'workers': 2, 'nprocs_per_endpoint': ranks, 'threads_per_rank': 1,
            'node_memory_GiB': memory_bytes/1024**3, 'allocated_cpus': cpus,
            'observed_prior_startup_max_rank_RSS_KiB': 32055664,
            'memory_allowance_GiB_per_rank': 64, 'node_memory_fraction_for_rank_planning': 0.75,
            'maximum_native_ranks_per_endpoint': 16}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cpus = _parse_allocation()
    memory = os.sysconf('SC_PHYS_PAGES')*os.sysconf('SC_PAGE_SIZE')
    m = read_json(args.manifest)
    if len(m['tasks']) != 2 or {t['metal'] for t in m['tasks']} != {'La','Ca'}:
        raise InvalidArtifact('this launch policy requires the approved paired endpoints')
    plan = layout(cpus, memory)
    plan.update(manifest=record(args.manifest), implementation=record(__file__),
                node=socket.gethostname(), job_id=os.environ['SLURM_JOB_ID'],
                compute_budget=None, wall_time_limit=None)
    write_new(args.output, plan)
    print(plan, flush=True)
    run_manifest(args.manifest, orca_path=verify(m['orca']), workers=2, nprocs=plan['nprocs_per_endpoint'])
