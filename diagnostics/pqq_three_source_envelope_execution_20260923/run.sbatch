#!/bin/bash
# Explicit finite plan;32CPU,oneH200,200000MiB;32x1 scalar workers.
set -uo pipefail
ENVELOPE_PLAN=$1
case "$ENVELOPE_PLAN" in /*) ;; *) exit 2;; esac
ENVELOPE_ROOT=$(dirname "$ENVELOPE_PLAN")
ENVELOPE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
ENVELOPE_SCRIPT="$ENVELOPE_ROOT/implementation/pqq_three_source_envelope_execution.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$ENVELOPE_PY" "$ENVELOPE_SCRIPT" execute --plan "$ENVELOPE_PLAN"
ENVELOPE_STATUS=$?
ENVELOPE_RESULT="$ENVELOPE_ROOT/RESULT_${SLURM_JOB_ID}.json"
"$ENVELOPE_PY" "$ENVELOPE_SCRIPT" collect --plan "$ENVELOPE_PLAN" --output "$ENVELOPE_RESULT"
ENVELOPE_COLLECT=$?
if [ "$ENVELOPE_COLLECT" -eq 0 ]; then
  "$ENVELOPE_PY" "$ENVELOPE_SCRIPT" report --result "$ENVELOPE_RESULT" --output "${ENVELOPE_RESULT%.json}.md"
  ENVELOPE_COLLECT=$?
fi
if [ "$ENVELOPE_STATUS" -ne 0 ]; then exit "$ENVELOPE_STATUS"; fi
exit "$ENVELOPE_COLLECT"
