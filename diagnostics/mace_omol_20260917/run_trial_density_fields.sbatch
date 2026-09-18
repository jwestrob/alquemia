#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-trial-density-fields
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_trial_density_fields.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TRIAL_FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
TRIAL_FIELD_ROOT=$(dirname "$1")
TRIAL_FIELD_DRIVER=$TRIAL_FIELD_ROOT/implementation/mace_trial_density_fields.py
"$TRIAL_FIELD_PY" "$TRIAL_FIELD_DRIVER" execute --manifest "$1" \
    --output "$TRIAL_FIELD_ROOT/execution_job_${SLURM_JOB_ID}.json"
"$TRIAL_FIELD_PY" "$TRIAL_FIELD_DRIVER" collect --manifest "$1" \
    --output "$TRIAL_FIELD_ROOT/report_job_${SLURM_JOB_ID}"
