#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-global-field
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: run_quantum.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FIELD_MANIFEST=$1
case "$FIELD_MANIFEST" in /*) ;; *) echo 'Manifest path must be absolute' >&2; exit 2;; esac
FIELD_OUTPUT_DIR=$(dirname "$FIELD_MANIFEST")
"$FIELD_PY" scripts/global_electrostatic.py execute --manifest "$FIELD_MANIFEST" --output "$FIELD_OUTPUT_DIR/execution_${SLURM_JOB_ID}.json"
"$FIELD_PY" scripts/global_electrostatic.py collect --manifest "$FIELD_MANIFEST" --output "$FIELD_OUTPUT_DIR/collection_${SLURM_JOB_ID}.json"
