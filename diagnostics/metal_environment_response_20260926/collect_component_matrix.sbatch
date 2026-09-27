#!/bin/bash
#SBATCH --job-name=Nikasha-MM-collect
#SBATCH --partition=standard-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --no-requeue
set -euo pipefail
FIELD_WORK="$1"
FIELD_JOB="$2"
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$FIELD_PY" "$FIELD_WORK/implementation/metal_environment_component_checks.py" collect --manifest "$FIELD_WORK/manifest.json" --output "$FIELD_WORK/FINAL_COLLECTION_${FIELD_JOB}.json" > "$FIELD_WORK/matrix_collection_${FIELD_JOB}_stdout.json"
sacct -j "$FIELD_JOB" --parsable2 --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,ReqMem,MaxRSS,ExitCode > "$FIELD_WORK/matrix_accounting_${FIELD_JOB}.txt"
