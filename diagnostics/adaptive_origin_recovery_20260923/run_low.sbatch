#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-origin-solvent
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
RECOVERY_POOL=$1
case "$RECOVERY_POOL" in /*) ;; *) exit 2;; esac
RECOVERY_RUN=$(dirname "$RECOVERY_POOL")
RECOVERY_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$RECOVERY_PY" "$RECOVERY_RUN/implementation/adaptive_origin_recovery.py" validate-pool --manifest "$RECOVERY_POOL"
RECOVERY_VALID=$?
if [ "$RECOVERY_VALID" -ne 0 ]; then exit "$RECOVERY_VALID"; fi
"$RECOVERY_PY" "$RECOVERY_RUN/implementation/affordable_workflow.py" execute --manifest "$RECOVERY_RUN/solvent/shard_0/manifest.json"
RECOVERY_EXEC=$?
"$RECOVERY_PY" "$RECOVERY_RUN/implementation/adaptive_origin_recovery.py" collect --manifest "$RECOVERY_POOL" --output "$RECOVERY_RUN/final_collection.json"
RECOVERY_COLLECT=$?
if [ "$RECOVERY_EXEC" -ne 0 ]; then exit "$RECOVERY_EXEC"; fi
exit "$RECOVERY_COLLECT"
