#!/bin/bash
#SBATCH -p gpu
#SBATCH -J precision-rank-2x2
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
RANK_READY=$1
case "$RANK_READY" in /*) ;; *) exit 2;; esac
RANK_ROOT=$(dirname "$RANK_READY")
RANK_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$RANK_PY" "$RANK_ROOT/implementation/affordable_workflow.py" execute --manifest "$RANK_ROOT/rank_1/manifest.json" &
RANK_ONE_PID=$!
"$RANK_PY" "$RANK_ROOT/implementation/affordable_workflow.py" execute --manifest "$RANK_ROOT/rank_8/manifest.json" &
RANK_EIGHT_PID=$!
wait "$RANK_ONE_PID"
RANK_ONE_STATUS=$?
wait "$RANK_EIGHT_PID"
RANK_EIGHT_STATUS=$?
"$RANK_PY" "$RANK_ROOT/implementation/precision_rank_geometry.py" collect --ready "$RANK_READY" --output "$RANK_ROOT/COLLECTION.json"
RANK_COLLECT_STATUS=$?
if [ "$RANK_ONE_STATUS" -ne 0 ]; then exit "$RANK_ONE_STATUS"; fi
if [ "$RANK_EIGHT_STATUS" -ne 0 ]; then exit "$RANK_EIGHT_STATUS"; fi
exit "$RANK_COLLECT_STATUS"
