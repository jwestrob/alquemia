#!/bin/bash
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --qos=standard_expand
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
exec /usr/bin/python3 "$1" "$2" "$3" "$4"
