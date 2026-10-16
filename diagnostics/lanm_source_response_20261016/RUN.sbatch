#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=LanM-saved-source-forces
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue

set -euo pipefail
: "${SLURM_JOB_ID:?Saved-force profiling requires Slurm}"
if [[ $# -ne 4 ]]; then
    echo 'Usage: RUN.sbatch ABS_PROFILE_PY ABS_PLAN_MD ABS_FINAL_PROVENANCE_JSON FRESH_ABS_OUTPUT_DIR' >&2
    exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  "$1" --plan "$2" --provenance "$3" --output "$4"
