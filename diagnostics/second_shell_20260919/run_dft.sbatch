#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J donor-second-shell-DFT
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
SECOND_MANIFEST=$1
case "$SECOND_MANIFEST" in /*) ;; *) exit 2;; esac
SECOND_OUT=$(dirname "$SECOND_MANIFEST")
SECOND_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
"$SECOND_PY" "$SECOND_OUT/implementation/affordable_workflow.py" execute --manifest "$SECOND_MANIFEST" --output "$SECOND_OUT/execution_${SLURM_JOB_ID}.json"
"$SECOND_PY" "$SECOND_OUT/implementation/second_shell_context.py" collect --manifest "$SECOND_MANIFEST" --output "$SECOND_OUT/collection_${SLURM_JOB_ID}.json"
