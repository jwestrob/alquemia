#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=LBT3-source-coordinates
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual allocation required}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python "$1" --repo-root "$2" --config "$3" --output "$4"
