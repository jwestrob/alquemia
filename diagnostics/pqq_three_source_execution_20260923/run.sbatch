#!/bin/bash
# HELD: use only after the separately recorded finite launch agreement.
# Existing profile: explicit H200 host,32CPUs,1GPU,200000MiB;32x1 scalar workers.
set -uo pipefail
THREE_PLAN=$1
case "$THREE_PLAN" in /*) ;; *) exit 2;; esac
THREE_ROOT=$(dirname "$THREE_PLAN")
THREE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
THREE_SCRIPT="$THREE_ROOT/implementation/pqq_three_source_execution.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$THREE_PY" "$THREE_SCRIPT" execute --plan "$THREE_PLAN"
THREE_STATUS=$?
THREE_RESULT="$THREE_ROOT/RESULT_${SLURM_JOB_ID}.json"
"$THREE_PY" "$THREE_SCRIPT" collect --plan "$THREE_PLAN" --output "$THREE_RESULT"
THREE_COLLECT_STATUS=$?
if [ "$THREE_COLLECT_STATUS" -eq 0 ]; then
    "$THREE_PY" "$THREE_SCRIPT" report --result "$THREE_RESULT" --output "${THREE_RESULT%.json}.md"
    THREE_COLLECT_STATUS=$?
fi
if [ "$THREE_STATUS" -ne 0 ]; then exit "$THREE_STATUS"; fi
exit "$THREE_COLLECT_STATUS"
