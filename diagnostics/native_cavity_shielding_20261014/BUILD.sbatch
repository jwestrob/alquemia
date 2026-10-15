#!/bin/bash
#SBATCH --job-name=native-cavity-build
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
exec /usr/bin/python3 "$(dirname "$1")/BUILD.py" "$1" "$2"
