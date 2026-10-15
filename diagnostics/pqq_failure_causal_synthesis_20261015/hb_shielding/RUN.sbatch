#!/bin/bash
#SBATCH --job-name=saved-HB-shielding
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Slurm required}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
HB_CONFIG="$1"
HB_SHA="$2"
HB_WORK=$(dirname "$HB_CONFIG")
HB_PYTHON=$(/usr/bin/python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["python"]["path"])' "$HB_CONFIG")
exec "$HB_PYTHON" "$HB_WORK/COMPARE.py" run --config "$HB_CONFIG" --config-sha "$HB_SHA"
