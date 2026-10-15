#!/bin/bash
#SBATCH --job-name=pqq-baseline-ordering
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Slurm required}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /usr/bin/python3 "$(dirname "$1")/REPORT.py" --config "$1" --config-sha "$2"
