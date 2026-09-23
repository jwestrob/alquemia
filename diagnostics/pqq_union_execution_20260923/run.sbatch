#!/bin/bash
# Submission remains pending root's integration/full-panel checkpoint.
# Both arms use the same explicit H200 host,32CPUs,1GPU,200000MiB.
set -uo pipefail
UNION_PLAN=$1
case "$UNION_PLAN" in /*) ;; *) exit 2;; esac
UNION_ROOT=$(dirname "$UNION_PLAN")
UNION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
UNION_SCRIPT="$UNION_ROOT/implementation/pqq_union_execution.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$UNION_PY" "$UNION_SCRIPT" execute --plan "$UNION_PLAN"
UNION_STATUS=$?
UNION_RESULT="$UNION_ROOT/RESULT.json"
if [ -e "$UNION_RESULT" ]; then UNION_RESULT="$UNION_ROOT/RESULT_${SLURM_JOB_ID}.json"; fi
"$UNION_PY" "$UNION_SCRIPT" collect --plan "$UNION_PLAN" --output "$UNION_RESULT"
UNION_COLLECT_STATUS=$?
if [ "$UNION_COLLECT_STATUS" -eq 0 ]; then
    "$UNION_PY" "$UNION_SCRIPT" report --result "$UNION_RESULT" --output "${UNION_RESULT%.json}.md"
    UNION_COLLECT_STATUS=$?
fi
if [ "$UNION_STATUS" -ne 0 ]; then exit "$UNION_STATUS"; fi
exit "$UNION_COLLECT_STATUS"
