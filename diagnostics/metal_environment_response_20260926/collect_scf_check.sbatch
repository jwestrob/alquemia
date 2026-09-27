#!/bin/bash
#SBATCH --job-name=Nikasha-scf-collect
#SBATCH --partition=standard-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --no-requeue
set -euo pipefail
FIELD_WORK="$1"
FIELD_JOB="$2"
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$FIELD_PY" "$FIELD_WORK/implementation/metal_environment_scf_check.py" collect \
    --manifest "$FIELD_WORK/manifest.json" --output "$FIELD_WORK/FINAL_COLLECTION.json" \
    > "$FIELD_WORK/final_collection_stdout.json"
sacct -j "$FIELD_JOB" --parsable2 --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,ReqMem,MaxRSS,ExitCode \
    > "$FIELD_WORK/final_accounting.txt"
"$FIELD_PY" - "$FIELD_WORK" "$FIELD_JOB" <<'PY'
import json,sys,datetime
from pathlib import Path
w=Path(sys.argv[1]); r=json.loads((w/'FINAL_COLLECTION.json').read_text())
counts={k:sum(c['status']==k for c in r['checks']) for k in ('passed','failed','unavailable')}
text=(f"# Native SCF diagnosis job {sys.argv[2]}: {r['status']}\n\n"
      f"Collected {datetime.datetime.now(datetime.timezone.utc).isoformat()}.\n\n"
      f"Checks: {counts}. Complete molecular cells: "
      f"{sum(v['status']=='complete' for v in r['rows'].values())}/{len(r['rows'])}.\n\n"
      "Inspect actual terminal SCF criteria and core/external rigid-transform "
      "checks in FINAL_COLLECTION.json. Qualification covers only the declared electronic "
      "component; new directional checks remain unrun. Full hybrid and ML comparison remain unsupported.\n")
with (w/'AUTO_REPORT.md').open('x') as f:f.write(text)
PY
