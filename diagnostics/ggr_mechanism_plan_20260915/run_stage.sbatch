#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-ggr
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [ "$#" -ne 2 ]; then
    echo 'Usage: run_stage.sbatch ABSOLUTE_MANIFEST energy|sensitivity' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
GGR_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
GGR_MANIFEST=$1
GGR_MODE=$2
GGR_OUTPUT_DIR=$(dirname "$GGR_MANIFEST")
case "$GGR_MANIFEST" in /*) ;; *) echo 'Manifest path must be absolute' >&2; exit 2;; esac
case "$GGR_MODE" in energy|sensitivity) ;; *) echo 'Unknown collector' >&2; exit 2;; esac
"$GGR_PY" scripts/affordable_workflow.py execute --manifest "$GGR_MANIFEST" --output "$GGR_OUTPUT_DIR/execution_${SLURM_JOB_ID}.json"
if [ "$GGR_MODE" = energy ]; then
    "$GGR_PY" scripts/ggr_workflow.py collect --manifest "$GGR_MANIFEST" --output "$GGR_OUTPUT_DIR/collection_${SLURM_JOB_ID}.json"
else
    "$GGR_PY" scripts/ggr_sensitivity.py collect --manifest "$GGR_MANIFEST" --output "$GGR_OUTPUT_DIR/collection_${SLURM_JOB_ID}.json"
fi
