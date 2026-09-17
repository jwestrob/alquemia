#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-vacuum-match
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_matched_vacuum.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
VACUUM_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
VACUUM_MANIFEST=$1
VACUUM_OUTPUT=$(dirname "$VACUUM_MANIFEST")
VACUUM_DRIVER="$VACUUM_OUTPUT/implementation/mace_omol_vacuum.py"
finish_collection() {
    VACUUM_EXIT=$?
    trap - EXIT
    "$VACUUM_PY" "$VACUUM_DRIVER" collect --manifest "$VACUUM_MANIFEST" \
        --output "$VACUUM_OUTPUT/collection_job_${SLURM_JOB_ID}.json" || true
    exit "$VACUUM_EXIT"
}
trap finish_collection EXIT
"$VACUUM_PY" "$VACUUM_DRIVER" execute --manifest "$VACUUM_MANIFEST" \
    --output "$VACUUM_OUTPUT/execution_${SLURM_JOB_ID}.json"
