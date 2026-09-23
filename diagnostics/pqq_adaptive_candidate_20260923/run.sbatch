#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-optin-adaptive
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
CANDIDATE_PLAN=$1
case "$CANDIDATE_PLAN" in /*) ;; *) exit 2;; esac
CANDIDATE_RUN=$(dirname "$CANDIDATE_PLAN")
CANDIDATE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
CANDIDATE_SCRIPT="$CANDIDATE_RUN/implementation/pqq_adaptive_candidate.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$CANDIDATE_PY" "$CANDIDATE_SCRIPT" execute --plan "$CANDIDATE_PLAN"
CANDIDATE_STATUS=$?
"$CANDIDATE_PY" "$CANDIDATE_SCRIPT" collect --plan "$CANDIDATE_PLAN" --output "$CANDIDATE_RUN/result.json"
CANDIDATE_COLLECT=$?
if [ "$CANDIDATE_COLLECT" -eq 0 ]; then
  "$CANDIDATE_PY" "$CANDIDATE_SCRIPT" report --result "$CANDIDATE_RUN/result.json" --output "$CANDIDATE_RUN/report.md"
fi
if [ "$CANDIDATE_STATUS" -ne 0 ]; then exit "$CANDIDATE_STATUS"; fi
exit "$CANDIDATE_COLLECT"
