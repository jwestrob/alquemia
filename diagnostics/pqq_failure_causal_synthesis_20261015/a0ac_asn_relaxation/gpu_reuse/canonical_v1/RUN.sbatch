#!/bin/bash
#SBATCH --job-name=canonical-MACE-GPU
#SBATCH --partition=gpu
#SBATCH --nodelist=node-128-512g-8gpu-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual Slurm required}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_THREAD_LIMIT=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
GPU_CANONICAL_CONFIG="$1"
GPU_CANONICAL_SHA="$2"
GPU_CANONICAL_WORK=$(dirname "$GPU_CANONICAL_CONFIG")
GPU_CANONICAL_PYTHON=$(/usr/bin/python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["python"]["path"])' "$GPU_CANONICAL_CONFIG")
exec "$GPU_CANONICAL_PYTHON" "$GPU_CANONICAL_WORK/ADAPTER.py" run --config "$GPU_CANONICAL_CONFIG" --config-sha "$GPU_CANONICAL_SHA"
