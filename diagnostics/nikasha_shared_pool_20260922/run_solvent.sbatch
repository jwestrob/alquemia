#!/bin/bash
#SBATCH -p standard
#SBATCH -J nikasha-pool-solvent
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
set -uo pipefail
POOL_MANIFEST=$1
POOL_SHARD=$2
case "$POOL_MANIFEST" in /*) ;; *) exit 2;; esac
case "$POOL_SHARD" in 0|1|2|3) ;; *) exit 2;; esac
POOL_RUN=$(dirname "$POOL_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
POOL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$POOL_PY" "$POOL_RUN/implementation/affordable_workflow.py" execute --manifest "$POOL_RUN/solvent/shard_$POOL_SHARD/manifest.json"
POOL_EXIT=$?
"$POOL_PY" "$POOL_RUN/implementation/nikasha_pool.py" collect --manifest "$POOL_MANIFEST" --output "$POOL_RUN/after_solvent_${POOL_SHARD}_${SLURM_JOB_ID}.json"
COLLECT_EXIT=$?
if [ "$POOL_EXIT" -ne 0 ]; then exit "$POOL_EXIT"; fi
exit "$COLLECT_EXIT"
