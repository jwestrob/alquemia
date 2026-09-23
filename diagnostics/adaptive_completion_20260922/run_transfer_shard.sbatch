#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-adaptive-fold
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
ADAPTIVE_MANIFEST=$1
ADAPTIVE_REFERENCE=$2
ADAPTIVE_SHARD=$3
case "$ADAPTIVE_MANIFEST" in /*) ;; *) exit 2;; esac
case "$ADAPTIVE_REFERENCE" in /*) ;; *) exit 2;; esac
case "$ADAPTIVE_SHARD" in 0|1|2|3) ;; *) exit 2;; esac
ADAPTIVE_RUN=$(dirname "$ADAPTIVE_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
ADAPTIVE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$ADAPTIVE_PY" "$ADAPTIVE_RUN/implementation/adaptive_completion_folds.py" execute --manifest "$ADAPTIVE_MANIFEST" --reference "$ADAPTIVE_REFERENCE" --shard-index "$ADAPTIVE_SHARD"
