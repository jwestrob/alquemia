#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-joint-proposals
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
JOINT_MANIFEST=$1
JOINT_POOL=$2
case "$JOINT_MANIFEST" in /*) ;; *) exit 2;; esac
case "$JOINT_POOL" in /*) ;; *) exit 2;; esac
JOINT_RUN=$(dirname "$JOINT_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
JOINT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$JOINT_PY" "$JOINT_RUN/implementation/adaptive_metal_proposals.py" execute --manifest "$JOINT_MANIFEST" --pool-collection "$JOINT_POOL"
"$JOINT_PY" "$JOINT_RUN/implementation/adaptive_metal_proposals.py" collect --manifest "$JOINT_MANIFEST" --output "$JOINT_RUN/after_proposals_${SLURM_JOB_ID}.json"
