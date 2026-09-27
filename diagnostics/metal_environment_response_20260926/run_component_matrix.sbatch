#!/bin/bash
#SBATCH --job-name=Nikasha-MM-checks
#SBATCH --partition=standard-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --no-requeue
set -euo pipefail
FIELD_WORK="$1"
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ "$(date +%s)" -ge "$(TZ=America/Los_Angeles date -d '2026-09-28 10:00:00' +%s)" ]; then
  exit 75
fi
"$FIELD_PY" "$FIELD_WORK/implementation/metal_environment_component_checks.py" execute --manifest "$FIELD_WORK/manifest.json" --workers 32 > "$FIELD_WORK/matrix_execution_${SLURM_JOB_ID}.json"
