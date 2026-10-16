#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=LBT3-saved-comparison
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Saved comparison requires Slurm}"
if [[ $# -ne 3 ]]; then
    echo 'Usage: RUN_COMPARE.sbatch ABS_COMPARE_PY ABS_CONFIG_JSON FRESH_ABS_OUTPUT_DIR' >&2
    exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python "$1" --config "$2" --output "$3"
