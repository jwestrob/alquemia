#!/bin/bash
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --qos=standard_expand
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual allocation required}"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
exec /usr/bin/python3 "$@"
