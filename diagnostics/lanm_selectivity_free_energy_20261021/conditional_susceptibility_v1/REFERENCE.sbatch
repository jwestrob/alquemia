#!/bin/bash
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --qos=standard_expand
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
set -euo pipefail
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
exec /home/jwestrob/.pyenv/versions/3.12.0/bin/python "$@"
