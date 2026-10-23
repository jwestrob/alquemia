#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard_expand
#SBATCH --account=banfield
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --no-requeue
#SBATCH --job-name=radial_field_recovery
set -euo pipefail
unset LD_PRELOAD LD_LIBRARY_PATH
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
exec /home/jwestrob/.pyenv/versions/3.12.0/bin/python3 "$1" run --config "$2" --out "$3"
