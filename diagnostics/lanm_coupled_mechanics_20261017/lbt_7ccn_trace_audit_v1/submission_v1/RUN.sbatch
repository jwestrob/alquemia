#!/bin/bash
#SBATCH --account=banfield
#SBATCH --qos=standard
#SBATCH --partition=gpu_h200
#SBATCH --job-name=lbt-saved-trace
#SBATCH --nodes=1
#SBATCH --nodelist=node-224-2t-8gpu-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual allocation required}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
exec "$@"
