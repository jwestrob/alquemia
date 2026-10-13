#!/bin/bash
#SBATCH --job-name=PedE-seed-check
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
if [ "$#" -ne 4 ]; then echo 'Usage: CHECK.sbatch dry-run|test|collect ABS_MANIFEST SHA ABS_OUTPUT' >&2;exit 2;fi
case "$1" in dry-run|test|collect) ;; *) exit 2;; esac
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
exec /usr/bin/python3 "$(dirname "$2")/implementation/SEEDED.py" "$@"
