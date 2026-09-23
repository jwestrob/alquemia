#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-union-collect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --nodelist=node-128-512g-8gpu-1
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNION_CPU_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$UNION_CPU_PY" - "$@" <<'PY'
import datetime, hashlib, json, os, pathlib, subprocess, sys, time
stage, reference, baseline, sources = [pathlib.Path(p).resolve() for p in sys.argv[1:]]
script = stage/'implementation/consistent_context_compare.py'
result = stage/'collection_final_v1.json'; comparison = stage/'comparison_v1.json'
commands = [
    [sys.executable, str(script), 'collect', '--stage', str(stage), '--output', str(result)],
    [sys.executable, str(script), 'transfer-compare', '--scope', 'primary225', '--collection', str(result),
     '--reference-path', str(reference), '--baseline', str(baseline), '--sources', str(sources), '--output', str(comparison)]]
started = time.monotonic(); receipts = []
for i, command in enumerate(commands):
    log = stage/f'autocollect_step_{i}.log'
    t = time.monotonic()
    with log.open('x') as f: run = subprocess.run(command, stdout=f, stderr=subprocess.STDOUT)
    receipts.append({'command': command, 'returncode': run.returncode, 'elapsed_seconds': time.monotonic()-t, 'log': str(log)})
    if run.returncode: break
record = {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'job_id': os.environ.get('SLURM_JOB_ID'), 'elapsed_seconds': time.monotonic()-started,
          'allocated_cpus': int(os.environ.get('SLURM_CPUS_PER_TASK', '1')), 'commands': receipts,
          'analysis_script': {'path': str(script), 'sha256': hashlib.sha256(script.read_bytes()).hexdigest()},
          'new_molecular_calls': 0, 'status': 'complete' if len(receipts) == 2 and not receipts[-1]['returncode'] else 'failed'}
with (stage/'autocollect_execution.json').open('x') as f: json.dump(record, f, indent=2); f.write('\n')
sys.exit(receipts[-1]['returncode'])
PY
