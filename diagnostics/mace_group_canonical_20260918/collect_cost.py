"""Collect actual job accounting and every molecular attempt, including failures."""
import argparse
import datetime
import json
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from mace_hybrid import accepted_attempt


def cpu_seconds(text):
    days, clock = text.split('-', 1) if '-' in text else ('0', text)
    value = 0.
    for part in clock.split(':'):
        value = value * 60 + float(part)
    return int(days) * 86400 + value


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('manifest', 'submission', 'preparation', 'output'):
        p.add_argument('--' + name, required=True)
    p.add_argument('--allow-running', action='store_true')
    p.add_argument('--preparation-reused', action='store_true',
                   help='Preparation timing is historical, not new work for this job.')
    p.add_argument('--prior-costs', default='Four cached crystal endpoints, source physical preparation and numerical/grouping qualification are additional historical costs.')
    a = p.parse_args(); mp = Path(a.manifest).resolve(); m = read_json(mp); submission = read_json(a.submission)
    if submission['manifest'] != record(mp) or submission['returncode'] != 0:
        raise InvalidArtifact('successful submission must match the exact manifest')
    job = str(submission['job_id'])
    result = subprocess.run(['sacct', '-j', job, '--parsable2', '--noheader',
                             '--format=JobIDRaw,State,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS,AllocTRES%100'],
                            text=True, capture_output=True, check=True)
    rows = [line.split('|') for line in result.stdout.splitlines() if line.strip()]
    parent = next(r for r in rows if r[0] == job)
    terminal = parent[1].split()[0] in {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL', 'PREEMPTED'}
    if not terminal and not a.allow_running:
        raise InvalidArtifact('job still live; final cost unavailable')
    args = submission['args']
    if '--gres=gpu:1' not in args:
        raise InvalidArtifact('this collector requires the declared one-GPU allocation')
    attempts = []; successful = {}; tasks = {t['task_id']: t for t in m['tasks']}
    for tid, t in tasks.items():
        for attempt in sorted((mp.parent / 'execution' / tid).glob('attempt_*')):
            path = attempt / 'result.json'; receipt = attempt / 'receipt.json'
            r = read_json(path) if path.exists() else {}
            accepted = accepted_attempt(attempt, t, mp) is not None
            if accepted:
                successful[tid] = r
            attempts.append({'task_id': tid, 'accepted': accepted, 'status': r.get('status', 'incomplete'),
                             'result': record(path) if path.exists() else None,
                             'receipt': record(receipt) if receipt.exists() else None,
                             'evaluation_seconds': r.get('evaluation_seconds'), 'worker_seconds': r.get('wall_seconds'),
                             'failure_reason': r.get('reason'), 'peak_cuda_allocated_bytes': r.get('peak_cuda_allocated_bytes'),
                             'peak_cuda_reserved_bytes': r.get('peak_cuda_reserved_bytes'),
                             'peak_worker_host_RSS_KiB': r.get('peak_host_RSS_KiB')})
    out = Path(a.output).resolve(); out.mkdir(parents=True, exist_ok=False)
    raw = out / 'sacct.tsv'; raw.write_text(result.stdout)
    prep = read_json(a.preparation)
    seconds = int(parent[2]); cores = int(parent[3])
    data = {'status': 'final' if terminal else 'partial', 'job': job, 'scheduler_state': parent[1],
            'collected_at_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'manifest': record(mp), 'submission': record(a.submission), 'sacct': record(raw),
            'successful_new_endpoints': len(successful), 'declared_new_endpoints': len(tasks),
            'failed_model_attempts': sum(r['status'] == 'failed' for r in attempts),
            'unfinished_attempts': sum(r['status'] == 'incomplete' for r in attempts),
            'elapsed_allocation_seconds': seconds, 'GPU_allocation_seconds': seconds,
            'allocated_core_seconds': seconds * cores, 'allocated_CPUs': cores,
            'reported_CPU_seconds': cpu_seconds(parent[4]) if terminal and parent[4] else None,
            'summed_model_evaluation_seconds': sum(r['evaluation_seconds'] or 0 for r in attempts),
            'summed_worker_seconds': sum(r['worker_seconds'] or 0 for r in attempts),
            'model_timing_missing_attempts': sum(r['evaluation_seconds'] is None for r in attempts),
            'preparation': record(a.preparation), 'preparation_wall_seconds': prep['wall_seconds'],
            'preparation_CPU_seconds': prep['CPU_seconds'], 'attempts': attempts,
            'preparation_reused': a.preparation_reused,
            'new_DFT_calls': 0, 'new_solvent_calls': 0,
            'prior_costs': a.prior_costs,
            'other_local_costs': 'Manifest preparation, failed/successful preflight, tests, reports and plotting are additional and not fully CPU-profiled.',
            'production_cost_established': False}
    for field in ('peak_cuda_allocated_bytes', 'peak_cuda_reserved_bytes', 'peak_worker_host_RSS_KiB'):
        data[field] = max((r[field] for r in attempts if r[field] is not None), default=None)
    source = out / 'collect_cost.py'; source.write_bytes(Path(__file__).read_bytes()); data['implementation'] = record(source)
    write_new(out / 'result.json', data)
    print(json.dumps({k: v for k, v in data.items() if k != 'attempts'}, indent=2))


if __name__ == '__main__':
    main()
