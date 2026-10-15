#!/bin/bash
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --job-name=ExaA-cavity
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual Slurm required}"
test "$#" -eq 2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec /usr/bin/python3 "$(dirname "$1")/RUN.py" execute --config "$1" --config-sha "$2"
