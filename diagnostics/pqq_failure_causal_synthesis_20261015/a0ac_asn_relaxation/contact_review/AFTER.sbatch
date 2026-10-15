#!/bin/bash
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=Asn-final-contacts
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual Slurm required}"
test "$#" -eq 7
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python - "$@" <<'PY'
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

source, collection, output, receipt, code_sha, method_sha, wrapper_sha = sys.argv[1:]
source, collection, output, receipt = map(Path, (source, collection, output, receipt))
if receipt.exists():
    raise FileExistsError(receipt)
record = {
    'job_id': os.environ['SLURM_JOB_ID'], 'status': 'failed',
    'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'collection': str(collection), 'output': str(output), 'new_energy_calls': 0,
    'scope': 'Terminal remaining-stage contact analysis only; no molecular execution',
}
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
try:
    expected = {'ANALYZE.py': code_sha, 'METHOD.md': method_sha, 'RUN.sbatch': wrapper_sha}
    for name, sha in expected.items():
        if digest(source / name) != sha:
            raise ValueError('Qualified source changed: ' + name)
    record['qualified_source_hashes'] = expected
    binding = subprocess.check_output([
        sys.executable, str(source / 'ANALYZE.py'), 'freeze',
        '--collection', str(collection), '--output', str(output),
    ], text=True)
    config = json.loads(binding)
    record['config'] = config
    if digest(Path(config['path'])) != config['sha256']:
        raise ValueError('Bound config changed')
    for name, sha in expected.items():
        if digest(output / name) != sha:
            raise ValueError('Copied source differs: ' + name)
    subprocess.run([
        sys.executable, str(output / 'ANALYZE.py'), 'execute',
        '--config', config['path'], '--config-sha', config['sha256'],
    ], check=True)
    execution = output / 'execution/EXECUTION.json'
    actual = json.loads(execution.read_text())
    record.update(status=actual['status'], execution={
        'path': str(execution), 'sha256': digest(execution)},
        all_candidate_geometries_available=actual['all_candidate_geometries_available'])
except BaseException:
    record['traceback'] = traceback.format_exc()
    raise
finally:
    record['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with receipt.open('x') as f:
        json.dump(record, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps(record))
PY
