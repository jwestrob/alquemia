#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-adaptive-angular
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
ADAPTIVE_MANIFEST=$1
ADAPTIVE_POOL=$2
case "$ADAPTIVE_MANIFEST" in /*) ;; *) exit 2;; esac
case "$ADAPTIVE_POOL" in /*) ;; *) exit 2;; esac
ADAPTIVE_RUN=$(dirname "$ADAPTIVE_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
ADAPTIVE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$ADAPTIVE_PY" "$ADAPTIVE_RUN/implementation/adaptive_angular_proposals.py" execute --manifest "$ADAPTIVE_MANIFEST" --pool-collection "$ADAPTIVE_POOL"
ADAPTIVE_EXIT=$?
"$ADAPTIVE_PY" "$ADAPTIVE_RUN/implementation/adaptive_angular_proposals.py" collect --manifest "$ADAPTIVE_MANIFEST" --output "$ADAPTIVE_RUN/after_proposals_${SLURM_JOB_ID}.json"
ADAPTIVE_COLLECT_EXIT=$?
if [ "$ADAPTIVE_EXIT" -ne 0 ]; then exit "$ADAPTIVE_EXIT"; fi
exit "$ADAPTIVE_COLLECT_EXIT"
