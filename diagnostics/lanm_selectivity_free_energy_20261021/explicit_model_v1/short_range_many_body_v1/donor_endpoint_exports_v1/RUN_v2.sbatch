#!/bin/bash
#SBATCH --job-name=srmb-donor8
#SBATCH --partition=standard
#SBATCH --qos=standard_expand
#SBATCH --account=banfield
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --no-requeue
set -euo pipefail
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
unset LD_PRELOAD LD_LIBRARY_PATH
exec /home/jwestrob/.pyenv/versions/3.12.0/bin/python "${1:?RUN.py}" --manifest "${2:?manifest}" --release "${3:?root release}" --out "${4:?fresh output}"
