#!/bin/bash
#SBATCH --job-name=PQQ-pose-prep
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
POSE_MANIFEST="$1"
POSE_SHA="$2"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
export OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
POSE_WORK=$(dirname "$POSE_MANIFEST")
POSE_PYTHON=$(/usr/bin/python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["python"]["path"])' "$POSE_MANIFEST")
exec "$POSE_PYTHON" "$POSE_WORK/implementation/PREPARE.py" execute --manifest "$POSE_MANIFEST" --manifest-sha "$POSE_SHA"
