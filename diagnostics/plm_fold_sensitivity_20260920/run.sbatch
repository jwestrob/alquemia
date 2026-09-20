#!/bin/bash
#SBATCH -p gpu
#SBATCH -J PLM-three-fold
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
PLM_STANDARD_PLAN=$1
case "$PLM_STANDARD_PLAN" in /*) ;; *) exit 2;; esac
PLM_STANDARD_RUN=$(dirname "$PLM_STANDARD_PLAN")
PLM_STANDARD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$PLM_STANDARD_PY" "$PLM_STANDARD_RUN/implementation/affordable_workflow.py" standard execute --plan "$PLM_STANDARD_PLAN"
PLM_EXECUTE_STATUS=$?
if [ "$PLM_EXECUTE_STATUS" -eq 0 ]; then
    "$PLM_STANDARD_PY" "$PLM_STANDARD_RUN/implementation/affordable_workflow.py" standard collect --plan "$PLM_STANDARD_PLAN" --output "$PLM_STANDARD_RUN/result_${SLURM_JOB_ID}.json"
    PLM_COLLECT_STATUS=$?
    if [ "$PLM_COLLECT_STATUS" -ne 0 ]; then exit "$PLM_COLLECT_STATUS"; fi
    "$PLM_STANDARD_PY" "$PLM_STANDARD_RUN/implementation/affordable_workflow.py" standard ensemble --plan "$PLM_STANDARD_PLAN" --collection "$PLM_STANDARD_RUN/result_${SLURM_JOB_ID}.json" --output "$PLM_STANDARD_RUN/ensemble_${SLURM_JOB_ID}.json"
    exit $?
else
    # Preserve current preparation failures and all required members; no retries.
    "$PLM_STANDARD_PY" "$PLM_STANDARD_RUN/implementation/affordable_workflow.py" standard ensemble --plan "$PLM_STANDARD_PLAN" --output "$PLM_STANDARD_RUN/ensemble_${SLURM_JOB_ID}.json"
    exit "$PLM_EXECUTE_STATUS"
fi
