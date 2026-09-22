#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-common-pool
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
POOL_MANIFEST=$1
case "$POOL_MANIFEST" in /*) ;; *) exit 2;; esac
POOL_RUN=$(dirname "$POOL_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
POOL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$POOL_PY" "$POOL_RUN/implementation/nikasha_pool.py" execute-mace --manifest "$POOL_MANIFEST"
"$POOL_PY" "$POOL_RUN/implementation/nikasha_pool.py" collect --manifest "$POOL_MANIFEST" --output "$POOL_RUN/after_MACE_${SLURM_JOB_ID}.json"
