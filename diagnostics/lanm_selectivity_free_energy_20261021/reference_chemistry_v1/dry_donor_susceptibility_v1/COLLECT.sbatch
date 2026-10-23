#!/bin/bash
#SBATCH --job-name=dry-field-collect
#SBATCH --qos=standard_expand
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
set -euo pipefail
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
exec /home/jwestrob/.pyenv/versions/3.12.0/bin/python "$1" collect --panel "$2" --inventory "$3" --stage "$4" --out "$5"
