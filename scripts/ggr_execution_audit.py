"""Collect actual GGR endpoint receipts and Slurm allocation costs; no execution."""
import argparse
from datetime import datetime
from pathlib import Path
import subprocess

from affordable_common import InvalidArtifact, read_json, record, verify, write_new


def audit(manifests, output):
    endpoints, jobs = [], set()
    for path in manifests:
        manifest = read_json(path)
        for task in manifest['tasks']:
            rp = Path(task['output_path'] + '.execution.json')
            item = {k: task.get(k) for k in ('task_id', 'task_type', 'case', 'representation', 'metal')}
            item.update(manifest=record(path), status='missing', wall_seconds=None,
                        assigned_rank_seconds=None, receipt=None)
            if rp.exists():
                receipt = read_json(rp)
                if receipt['manifest'] != record(path) or receipt['task_id'] != task['task_id']:
                    raise InvalidArtifact('receipt task/manifest mismatch')
                for artifact in receipt['artifacts'].values():
                    verify(artifact)
                wall = (datetime.fromisoformat(receipt['finished_at_utc']) -
                        datetime.fromisoformat(receipt['started_at_utc'])).total_seconds()
                if wall < 0:
                    raise InvalidArtifact('negative receipt duration')
                allocation = receipt['allocation']
                jid = allocation['slurm_job_id']
                if not str(jid).isdigit():
                    raise InvalidArtifact('missing/unexpected Slurm job identity')
                jobs.add(str(jid))
                complete = receipt['normal_termination'] and receipt['scf_converged'] and receipt['returncode'] == 0
                ranks = receipt['parallelism']['nprocs']
                item.update(receipt=record(rp), status='complete' if complete else 'failed',
                            wall_seconds=wall, mpi_ranks=ranks, assigned_rank_seconds=wall*ranks,
                            allocation=allocation, orca_executable=receipt['orca_executable'])
            endpoints.append(item)
    command = ['sacct', '-j', ','.join(sorted(jobs)), '--format',
               'JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,TotalCPU,MaxRSS', '-P']
    if not jobs:
        raise InvalidArtifact('no actual Slurm receipts available')
    result = subprocess.run(command, text=True, capture_output=True, check=True)
    lines = result.stdout.strip().splitlines()
    header = lines[0].split('|')
    accounting = [dict(zip(header, line.split('|'))) for line in lines[1:]]
    parent = [r for r in accounting if r['JobID'] in jobs]
    if len(parent) != len(jobs):
        raise InvalidArtifact('incomplete Slurm parent accounting')
    payload = {'schema_version': 'alquemia.ggr_execution_audit.v1', 'implementation': record(__file__),
               'manifests': [record(p) for p in manifests], 'endpoints': endpoints,
               'completed_endpoints': sum(e['status'] == 'complete' for e in endpoints),
               'failed_endpoints': sum(e['status'] == 'failed' for e in endpoints),
               'missing_endpoints': sum(e['status'] == 'missing' for e in endpoints),
               'total_endpoints': len(endpoints),
               'sum_endpoint_wall_seconds': sum(e['wall_seconds'] or 0 for e in endpoints),
               'sum_endpoint_assigned_rank_seconds': sum(e['assigned_rank_seconds'] or 0 for e in endpoints),
               'allocated_core_seconds': sum(int(r['CPUTimeRAW']) for r in parent),
               'all_jobs_complete': all(r['State'] == 'COMPLETED' for r in parent),
               'sacct_command': command, 'sacct_raw': result.stdout, 'slurm_records': accounting,
               'gpu_seconds': 0, 'compute_budget': None, 'wall_time_limit': None,
               'limits': 'Endpoint sums count only available receipts, with missing and failed denominators separate. Slurm allocation cost includes idle ranks and wrapper overhead. Preparation/software-test/source-download costs are outside these batch receipts; unmeasured costs are unavailable.'}
    write_new(output, payload)
    return payload


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.manifest, args.output)
