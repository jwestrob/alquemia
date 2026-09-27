#!/bin/bash
#SBATCH --job-name=Nikasha-field-collect
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
"$FIELD_PY" "$FIELD_WORK/implementation/metal_environment_reference.py" collect \
    --manifest "$FIELD_WORK/manifest.json" --output "$FIELD_WORK/FINAL_COLLECTION_1219207.json" \
    > "$FIELD_WORK/final_collection_stdout_1219207.json"
sacct -j "$FIELD_JOB" --parsable2 --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,ReqMem,MaxRSS,ExitCode \
    > "$FIELD_WORK/final_accounting_1219207.txt"
# Independent of this chat/login process. Missing/failed cells stay missing.
"$FIELD_PY" - "$FIELD_WORK" "$FIELD_JOB" <<'PY'
import json,sys,datetime
from pathlib import Path
w=Path(sys.argv[1]); r=json.loads((w/'FINAL_COLLECTION_1219207.json').read_text())
rows=[f"| {k} | {v['status']} | {v.get('reason','')} |" for k,v in r['rows'].items()]
text=(f"# Native reference scout job {sys.argv[2]}: {r['status']}\n\n"
      f"Collected {datetime.datetime.now(datetime.timezone.utc).isoformat()}.\n\n"
      "| Cell | Status | Reason |\n|---|---|---|\n"+'\n'.join(rows)+'\n\n'
      f"Electronic double difference (kcal/mol): {r['delta_env_el_kcal_mol']}.\n\n"
      "This is an electronic embedding component diagnostic, not a complete QM/MM score. "
      "Native derivative checks, expanded-region qualification and ML comparison remain pending. "
      "No automatic scientific expansion or production promotion.\n\n"
      f"Actual results/accounting: `{w}/FINAL_COLLECTION_1219207.json`, `{w}/final_accounting_1219207.txt`.\n")
with (w/'AUTO_REPORT_1219207.md').open('x') as f:f.write(text)
PY
