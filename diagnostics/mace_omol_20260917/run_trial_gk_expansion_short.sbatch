#!/bin/bash
#SBATCH -p gpu
#SBATCH -J alquemia-trial-GK-short
#SBATCH --nodes=1
#SBATCH --nodelist=node-128-512g-8gpu-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=64474M
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_trial_gk_expansion_short.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
TRIAL_SHORT_MANIFEST=$1
TRIAL_SHORT_ROOT=$(dirname "$TRIAL_SHORT_MANIFEST")
TRIAL_SHORT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
TRIAL_SHORT_DRIVER="$TRIAL_SHORT_ROOT/implementation/mace_hybrid.py"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK" MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export PYTHONUNBUFFERED=1
finish_collection() {
    TRIAL_SHORT_EXIT=$?
    trap - EXIT
    "$TRIAL_SHORT_PY" "$TRIAL_SHORT_DRIVER" collect --manifest "$TRIAL_SHORT_MANIFEST" \
        --output "$TRIAL_SHORT_ROOT/collection_job_${SLURM_JOB_ID}.json" || true
    exit "$TRIAL_SHORT_EXIT"
}
trap finish_collection EXIT
"$TRIAL_SHORT_PY" "$TRIAL_SHORT_DRIVER" execute --manifest "$TRIAL_SHORT_MANIFEST" --memory-mode native
