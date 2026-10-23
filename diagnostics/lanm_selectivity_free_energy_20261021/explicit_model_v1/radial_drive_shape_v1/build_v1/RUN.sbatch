#!/bin/bash
#SBATCH --partition=standard_expand
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=radial_shape_build
set -euo pipefail
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
exec /home/jwestrob/.pyenv/versions/3.12.0/bin/python3 "$1" "$2" "$3"
