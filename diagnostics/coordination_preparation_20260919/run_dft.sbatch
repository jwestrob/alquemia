#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J coordination-context-DFT
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
COORD_MANIFEST=$1
case "$COORD_MANIFEST" in /*) ;; *) exit 2;; esac
COORD_OUT=$(dirname "$COORD_MANIFEST")
COORD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
"$COORD_PY" "$COORD_OUT/implementation/affordable_workflow.py" execute --manifest "$COORD_MANIFEST" --output "$COORD_OUT/execution_${SLURM_JOB_ID}.json"
"$COORD_PY" "$COORD_OUT/implementation/coordination_preparation_context.py" collect --manifest "$COORD_MANIFEST" --output "$COORD_OUT/collection_${SLURM_JOB_ID}.json"
