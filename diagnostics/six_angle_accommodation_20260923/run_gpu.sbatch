#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-six-angle
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
SIX_MANIFEST=$1
case "$SIX_MANIFEST" in /*) ;; *) exit 2;; esac
SIX_RUN=$(dirname "$SIX_MANIFEST")
SIX_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
SIX_SCRIPT="$SIX_RUN/implementation/six_angle_accommodation.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$SIX_PY" "$SIX_SCRIPT" execute --manifest "$SIX_MANIFEST"
"$SIX_PY" "$SIX_SCRIPT" prepare_low --manifest "$SIX_MANIFEST" --gpu-collection "$SIX_RUN/GPU_COLLECTION.json" --output "$SIX_RUN/scalar"
