#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=LBT-local-GFN
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual allocation required}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python "$1" --config "$2" --output "$3"
