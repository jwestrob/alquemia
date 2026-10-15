#!/bin/bash
#SBATCH --job-name=amide-transfer-MACE
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
TRANSFER_MACE_CONFIG="$1"
TRANSFER_MACE_SHA="$2"
TRANSFER_MACE_ACTION="${3:-run}"
TRANSFER_MACE_WORK=$(dirname "$TRANSFER_MACE_CONFIG")
TRANSFER_MACE_PYTHON=$(/usr/bin/python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["python"]["path"])' "$TRANSFER_MACE_CONFIG")
exec "$TRANSFER_MACE_PYTHON" "$TRANSFER_MACE_WORK/MACE.py" "$TRANSFER_MACE_ACTION" --config "$TRANSFER_MACE_CONFIG" --config-sha "$TRANSFER_MACE_SHA"
