#!/bin/bash
#SBATCH --job-name=MACE-GPU-reuse
#SBATCH --partition=gpu,gpu_h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual Slurm required}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_THREAD_LIMIT=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
GPU_SCOUT_CONFIG="$1"
GPU_SCOUT_SHA="$2"
GPU_SCOUT_WORK=$(dirname "$GPU_SCOUT_CONFIG")
GPU_SCOUT_PYTHON=$(/usr/bin/python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["parent_python"]["path"])' "$GPU_SCOUT_CONFIG")
exec "$GPU_SCOUT_PYTHON" "$GPU_SCOUT_WORK/SCOUT.py" run --config "$GPU_SCOUT_CONFIG" --config-sha "$GPU_SCOUT_SHA"
