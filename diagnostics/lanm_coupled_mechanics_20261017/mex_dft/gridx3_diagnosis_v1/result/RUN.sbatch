#!/bin/bash
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual allocation required}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11 "$@"
