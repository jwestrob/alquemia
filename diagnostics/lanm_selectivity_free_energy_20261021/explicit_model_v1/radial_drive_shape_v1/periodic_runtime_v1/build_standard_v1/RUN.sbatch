#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=banfield
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --no-requeue
#SBATCH --job-name=periodic_v2_build
set -euo pipefail
unset LD_PRELOAD LD_LIBRARY_PATH FFLAGS LDFLAGS CMAKE_GENERATOR CMAKE_TOOLCHAIN_FILE
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
exec /home/jwestrob/.pyenv/versions/3.12.0/bin/python3 "$1" --config "$2" --release "$3" --out "$4"
