#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J water-DFT-adjudication
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: run.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
HYDRATION_MANIFEST=$1
case "$HYDRATION_MANIFEST" in /*) ;; *) echo 'Manifest must be absolute' >&2; exit 2;; esac
HYDRATION_OUTPUT_DIR=$(dirname "$HYDRATION_MANIFEST")
HYDRATION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
"$HYDRATION_PY" scripts/affordable_workflow.py execute --manifest "$HYDRATION_MANIFEST" --output "$HYDRATION_OUTPUT_DIR/execution_${SLURM_JOB_ID}.json"
PYTHONPATH="$HYDRATION_OUTPUT_DIR/implementation:$HYDRATION_OUTPUT_DIR/implementation_${SLURM_JOB_ID}" "$HYDRATION_PY" "$HYDRATION_OUTPUT_DIR/implementation/hydration_adjudicate.py" collect --manifest "$HYDRATION_MANIFEST" --output "$HYDRATION_OUTPUT_DIR/collection_${SLURM_JOB_ID}.json"
