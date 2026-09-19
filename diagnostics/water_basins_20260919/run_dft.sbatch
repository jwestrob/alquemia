#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J water-finite-native
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: run_dft.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
WATER_BASIN_MANIFEST=$1
case "$WATER_BASIN_MANIFEST" in /*) ;; *) echo 'Manifest must be absolute' >&2; exit 2;; esac
WATER_BASIN_OUTPUT=$(dirname "$WATER_BASIN_MANIFEST")
WATER_BASIN_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
"$WATER_BASIN_PY" scripts/affordable_workflow.py execute --manifest "$WATER_BASIN_MANIFEST" --output "$WATER_BASIN_OUTPUT/execution_${SLURM_JOB_ID}.json"
PYTHONPATH="$WATER_BASIN_OUTPUT/implementation:$WATER_BASIN_OUTPUT/implementation_${SLURM_JOB_ID}" "$WATER_BASIN_PY" "$WATER_BASIN_OUTPUT/implementation/water_basin_validation.py" collect --manifest "$WATER_BASIN_MANIFEST" --output "$WATER_BASIN_OUTPUT/collection_${SLURM_JOB_ID}.json"
