#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-trial-GK-structures
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_trial_gk_expansion_quantum.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
TRIAL_GK_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
TRIAL_GK_MANIFEST=$1
TRIAL_GK_OUTPUT=$(dirname "$TRIAL_GK_MANIFEST")
TRIAL_GK_DRIVER="$TRIAL_GK_OUTPUT/implementation/mace_density_inputs.py"
finish_collection() {
    TRIAL_GK_EXIT=$?
    trap - EXIT
    "$TRIAL_GK_PY" "$TRIAL_GK_DRIVER" collect --manifest "$TRIAL_GK_MANIFEST" \
        --output "$TRIAL_GK_OUTPUT/collection_job_${SLURM_JOB_ID}.json" || true
    exit "$TRIAL_GK_EXIT"
}
trap finish_collection EXIT
"$TRIAL_GK_PY" "$TRIAL_GK_DRIVER" execute --manifest "$TRIAL_GK_MANIFEST" \
    --output "$TRIAL_GK_OUTPUT/execution_${SLURM_JOB_ID}.json"
