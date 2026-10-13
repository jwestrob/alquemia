#!/bin/bash
#SBATCH --job-name=PedE-DFT-parser-recovery
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
if [ "$#" -ne 2 ]; then echo 'Usage: RUN.sbatch ABS_CONFIG SHA' >&2; exit 2; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
exec /usr/bin/python3 "$(dirname "$1")/implementation/RECOVER.py" "$1" --sha "$2"
