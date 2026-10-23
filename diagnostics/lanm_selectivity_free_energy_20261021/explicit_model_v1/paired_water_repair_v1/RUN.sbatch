#!/bin/bash
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --qos=standard_expand
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual allocation required}"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
unset OMP_PLACES OMP_PROC_BIND
exec /usr/bin/python3 "$@"
