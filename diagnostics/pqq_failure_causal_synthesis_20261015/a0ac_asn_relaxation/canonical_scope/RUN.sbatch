#!/bin/bash
#SBATCH --job-name=canonical-amide-scope
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_THREAD_LIMIT=1
exec /usr/bin/python3 "$(dirname "$1")/INVENTORY.py" execute --manifest "$1" --sha "$2"
