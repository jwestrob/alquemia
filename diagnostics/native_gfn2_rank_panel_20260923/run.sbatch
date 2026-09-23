#!/bin/bash
# Explicit32CPU/64GiB CPU-only allocation is pinned by the submission receipt.
set -uo pipefail
PANEL_MANIFEST=$1
case "$PANEL_MANIFEST" in /*) ;; *) exit 2;; esac
PANEL_ROOT=$(dirname "$PANEL_MANIFEST")
PANEL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
PANEL_SCRIPT="$PANEL_ROOT/implementation/native_gfn2_rank_panel.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$PANEL_PY" "$PANEL_SCRIPT" execute --manifest "$PANEL_MANIFEST"
PANEL_STATUS=$?
"$PANEL_PY" "$PANEL_SCRIPT" collect --manifest "$PANEL_MANIFEST" --output "$PANEL_ROOT/collection.json"
PANEL_COLLECT_STATUS=$?
if [ "$PANEL_COLLECT_STATUS" -eq 0 ]; then
    "$PANEL_PY" "$PANEL_SCRIPT" compare --collection "$PANEL_ROOT/collection.json" --output "$PANEL_ROOT/COMPARISON.json"
    PANEL_COLLECT_STATUS=$?
fi
if [ "$PANEL_STATUS" -ne 0 ]; then exit "$PANEL_STATUS"; fi
exit "$PANEL_COLLECT_STATUS"
