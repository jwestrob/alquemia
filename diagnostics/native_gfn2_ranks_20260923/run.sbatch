#!/bin/bash
# Resources and afterany sequencing are declared in the saved submission receipt.
set -uo pipefail
RANK_MANIFEST=$1
RANK_DESIGN=$2
case "$RANK_MANIFEST" in /*) ;; *) exit 2;; esac
RANK_ROOT=$(dirname "$RANK_DESIGN")
RANK_DIRECTORY=$(dirname "$RANK_MANIFEST")
RANK_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
RANK_SCRIPT="$RANK_ROOT/implementation/native_gfn2_rank_scaling.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$RANK_PY" "$RANK_SCRIPT" execute --manifest "$RANK_MANIFEST"
RANK_STATUS=$?
"$RANK_PY" "$RANK_SCRIPT" collect --manifest "$RANK_MANIFEST" --output "$RANK_DIRECTORY/collection.json"
RANK_COLLECT_STATUS=$?
if [ "$RANK_DIRECTORY" = "$RANK_ROOT/rank_8" ]; then
 "$RANK_PY" "$RANK_SCRIPT" compare --design "$RANK_DESIGN" --output "$RANK_ROOT/COMPARISON.json"
fi
if [ "$RANK_STATUS" -ne 0 ]; then exit "$RANK_STATUS"; fi
exit "$RANK_COLLECT_STATUS"
