#!/bin/bash
#SBATCH --job-name=dry-component16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=standard_expand
set -euo pipefail
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE PYTHONDONTWRITEBYTECODE=1
exec /home/jwestrob/.pyenv/versions/3.12.0/bin/python "$1" run --manifest "$2" --release "$3" --out "$4"
