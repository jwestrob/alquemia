#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J response-core-DFT
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
RESPONSE_MANIFEST=$1
case "$RESPONSE_MANIFEST" in /*) ;; *) exit 2;; esac
RESPONSE_DIR=$(dirname "$RESPONSE_MANIFEST")
RESPONSE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
"$RESPONSE_PY" scripts/affordable_workflow.py execute --manifest "$RESPONSE_MANIFEST" --output "$RESPONSE_DIR/execution_${SLURM_JOB_ID}.json"
