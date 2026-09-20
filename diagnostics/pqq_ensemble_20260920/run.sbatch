#!/bin/bash
#SBATCH -p gpu
#SBATCH -J PQQ-standard
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
PQQ_STANDARD_PLAN=$1
case "$PQQ_STANDARD_PLAN" in /*) ;; *) exit 2;; esac
PQQ_STANDARD_RUN=$(dirname "$PQQ_STANDARD_PLAN")
PQQ_STANDARD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$PQQ_STANDARD_PY" "$PQQ_STANDARD_RUN/implementation/affordable_workflow.py" standard execute --plan "$PQQ_STANDARD_PLAN"
"$PQQ_STANDARD_PY" "$PQQ_STANDARD_RUN/implementation/affordable_workflow.py" standard collect --plan "$PQQ_STANDARD_PLAN" --output "$PQQ_STANDARD_RUN/result_${SLURM_JOB_ID}.json"
"$PQQ_STANDARD_PY" "$PQQ_STANDARD_RUN/implementation/affordable_workflow.py" standard ensemble --plan "$PQQ_STANDARD_PLAN" --collection "$PQQ_STANDARD_RUN/result_${SLURM_JOB_ID}.json" --output "$PQQ_STANDARD_RUN/ensemble_${SLURM_JOB_ID}.json"
