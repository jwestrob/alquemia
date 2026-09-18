#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-trial-GK-density
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_trial_gk_expansion_observations.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
TRIAL_OBS_MANIFEST=$1
TRIAL_OBS_ROOT=$(dirname "$TRIAL_OBS_MANIFEST")
TRIAL_OBS_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
TRIAL_OBS_DRIVER="$TRIAL_OBS_ROOT/implementation/mace_density_observations.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1
finish_collection() {
    TRIAL_OBS_EXIT=$?
    trap - EXIT
    "$TRIAL_OBS_PY" "$TRIAL_OBS_DRIVER" collect --manifest "$TRIAL_OBS_MANIFEST" \
        --output "$TRIAL_OBS_ROOT/report_job_${SLURM_JOB_ID}" || true
    exit "$TRIAL_OBS_EXIT"
}
trap finish_collection EXIT
"$TRIAL_OBS_PY" "$TRIAL_OBS_DRIVER" execute --manifest "$TRIAL_OBS_MANIFEST" \
    --output "$TRIAL_OBS_ROOT/execution_job_${SLURM_JOB_ID}.json"
