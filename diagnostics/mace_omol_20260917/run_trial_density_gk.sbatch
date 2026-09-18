#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-trial-density-gk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64000M
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_trial_density_gk.sbatch ABSOLUTE_PIPELINE_CONFIG' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=64
TRIAL_GK_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
TRIAL_GK_ROOT=$(dirname "$1")
TRIAL_GK_DRIVER=$TRIAL_GK_ROOT/implementation/mace_density_gk_hybrid.py
"$TRIAL_GK_PY" "$TRIAL_GK_DRIVER" pipeline-dry-run --config "$1"
"$TRIAL_GK_PY" "$TRIAL_GK_DRIVER" run-pipeline --config "$1"
