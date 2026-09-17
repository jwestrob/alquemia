#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-mechanics
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_quantum.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
MECHANICS_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
MECHANICS_MANIFEST=$1
MECHANICS_OUTPUT=$(dirname "$MECHANICS_MANIFEST")
finish_collection() {
    MECHANICS_EXIT=$?
    trap - EXIT
    "$MECHANICS_PY" scripts/mace_mechanics.py collect-dft --manifest "$MECHANICS_MANIFEST" \
        --output "$MECHANICS_OUTPUT/collection_job_${SLURM_JOB_ID}.json" || true
    exit "$MECHANICS_EXIT"
}
trap finish_collection EXIT
"$MECHANICS_PY" scripts/affordable_workflow.py execute --manifest "$MECHANICS_MANIFEST" \
    --output "$MECHANICS_OUTPUT/execution_${SLURM_JOB_ID}.json"
