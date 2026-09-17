#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-responsive-field
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_responsive_quantum.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
RESPONSIVE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
RESPONSIVE_MANIFEST=$1
RESPONSIVE_OUTPUT=$(dirname "$RESPONSIVE_MANIFEST")
RESPONSIVE_DRIVER="$RESPONSIVE_OUTPUT/implementation/mace_responsive_field.py"
finish_collection() {
    RESPONSIVE_EXIT=$?
    trap - EXIT
    "$RESPONSIVE_PY" "$RESPONSIVE_DRIVER" collect --manifest "$RESPONSIVE_MANIFEST" \
        --output "$RESPONSIVE_OUTPUT/collection_job_${SLURM_JOB_ID}.json" || true
    exit "$RESPONSIVE_EXIT"
}
trap finish_collection EXIT
"$RESPONSIVE_PY" "$RESPONSIVE_DRIVER" execute --manifest "$RESPONSIVE_MANIFEST" \
    --output "$RESPONSIVE_OUTPUT/execution_${SLURM_JOB_ID}.json"
