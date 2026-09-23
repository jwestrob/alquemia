#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-q88-repeat
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
REPEAT_MANIFEST=$1
case "$REPEAT_MANIFEST" in /*) ;; *) exit 2;; esac
REPEAT_RUN=$(dirname "$REPEAT_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
REPEAT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$REPEAT_PY" - "$REPEAT_MANIFEST" <<'PY'
import json,sys
from pathlib import Path
m=json.load(open(sys.argv[1]));dirs=[Path(t['output_path']).parent for t in m['tasks']]
assert len(dirs)==8 and all(not list(d.glob('*.gbw')) and not list(d.glob('*.xtbw')) for d in dirs), 'Fresh no-restart directories required'
print('All eight fresh directories have no GBW or xtbw restart artifacts',flush=True)
PY
REPEAT_CHECK=$?
if [ "$REPEAT_CHECK" -ne 0 ]; then exit "$REPEAT_CHECK"; fi
"$REPEAT_PY" "$REPEAT_RUN/implementation/affordable_workflow.py" execute --manifest "$REPEAT_MANIFEST"
REPEAT_EXIT=$?
"$REPEAT_PY" "$REPEAT_RUN/implementation/structure_informed_starts.py" collect-numeric --manifest "$REPEAT_MANIFEST" --output "$REPEAT_RUN/collection_${SLURM_JOB_ID}.json"
REPEAT_COLLECT_EXIT=$?
if [ "$REPEAT_EXIT" -ne 0 ]; then exit "$REPEAT_EXIT"; fi
exit "$REPEAT_COLLECT_EXIT"
