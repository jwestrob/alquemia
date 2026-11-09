#!/bin/bash
# Submit only after the declared periodic/throughput checks have been reviewed.
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JAX_PLATFORMS=cuda,cpu JAX_ENABLE_X64=true
export XLA_PYTHON_CLIENT_PREALLOCATE=false
python - <<'PY'
from pathlib import Path
import json,hashlib,os,datetime
assert os.environ.get('SLURM_JOB_ID')
d=Path('diagnostics/so3lr_neutral_pair_20261107')
w=Path('workspaces/so3lr_neutral_pair_20261107')
release=w/'MD_SCOUT_RELEASE.json'
assert release.exists(), 'Root review of periodic controls and throughput is required.'
r=json.loads(release.read_text());assert r['release'] is True
out=w/'md_scout_v1';out.mkdir(exist_ok=False)
files=[d/'PLAN.md',d/'MD_SCOUT.md',d/'HYDRATION_PREPARATION.md',d/'md_scout.yaml',d/'md_scout.sh',release,w/'hydrated_preparation_v1/RESULT.json',w/'hydrated_preparation_v1/start0/A.extxyz']
pins=[]
for p in files:
 b=p.read_bytes();(out/p.name).write_bytes(b)
 pins.append(dict(path=str(p.resolve()),sha256=hashlib.sha256(b).hexdigest()))
(out/'START.json').write_text(json.dumps(dict(job=os.environ['SLURM_JOB_ID'],created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),inputs=pins),indent=2)+'\n')
PY
timeout --signal=TERM --kill-after=10s 3600s srun -n1 -c12 --cpu-bind=cores \
 workspaces/so3lr_tauri_screen_20261106/software/env/bin/python -u -m so3lr.cli.so3lr_cli \
 --settings diagnostics/so3lr_neutral_pair_20261107/md_scout.yaml
