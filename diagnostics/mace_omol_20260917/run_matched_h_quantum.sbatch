#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-matched-H
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_matched_h_quantum.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
MATCHED_H_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
MATCHED_H_MANIFEST=$1
MATCHED_H_OUTPUT=$(dirname "$MATCHED_H_MANIFEST")
MATCHED_H_DRIVER="$MATCHED_H_OUTPUT/implementation/mace_omol_matched_h.py"
finish_collection() {
    MATCHED_H_EXIT=$?
    trap - EXIT
    "$MATCHED_H_PY" "$MATCHED_H_DRIVER" collect-quantum --manifest "$MATCHED_H_MANIFEST" \
        --output "$MATCHED_H_OUTPUT/collection_job_${SLURM_JOB_ID}.json" || true
    exit "$MATCHED_H_EXIT"
}
trap finish_collection EXIT
"$MATCHED_H_PY" "$MATCHED_H_DRIVER" execute-quantum --manifest "$MATCHED_H_MANIFEST" \
    --output "$MATCHED_H_OUTPUT/execution_${SLURM_JOB_ID}.json"
