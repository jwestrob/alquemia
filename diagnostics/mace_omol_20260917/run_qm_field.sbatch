#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-qm-electric-field
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_qm_field.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
QM_FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
QM_FIELD_ROOT=$(dirname "$1")
QM_FIELD_DRIVER=$QM_FIELD_ROOT/implementation/mace_qm_field.py
"$QM_FIELD_PY" "$QM_FIELD_DRIVER" execute --manifest "$1" \
    --output "$QM_FIELD_ROOT/execution_job_${SLURM_JOB_ID}.json"
"$QM_FIELD_PY" "$QM_FIELD_DRIVER" collect --manifest "$1" \
    --output "$QM_FIELD_ROOT/report_job_${SLURM_JOB_ID}"
