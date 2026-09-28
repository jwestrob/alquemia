#!/bin/bash
#SBATCH --job-name=Nikasha-Dy-collect
#SBATCH --partition=standard-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
set -euo pipefail
FIELD_WORK="$1"
FIELD_JOB="$2"
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$FIELD_PY" "$FIELD_WORK/implementation/metal_environment_lady_compact.py" collect --manifest "$FIELD_WORK/manifest.json" --output "$FIELD_WORK/FINAL_COLLECTION.json" > "$FIELD_WORK/collection_stdout.json"
sacct -j "$FIELD_JOB" --parsable2 --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,ReqMem,MaxRSS,ExitCode > "$FIELD_WORK/final_accounting.txt"
