#!/bin/bash
#SBATCH --job-name=PQQ-pose-MACE
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=7
#SBATCH --no-requeue
set -euo pipefail
POSE_CONFIG="$1"
POSE_SHA="$2"
POSE_ACTION="${3:-run}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
export OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
POSE_WORK=$(dirname "$POSE_CONFIG")
POSE_PYTHON=$(/usr/bin/python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["python"]["path"])' "$POSE_CONFIG")
exec "$POSE_PYTHON" "$POSE_WORK/SCORE.py" "$POSE_ACTION" --config "$POSE_CONFIG" --config-sha "$POSE_SHA"
