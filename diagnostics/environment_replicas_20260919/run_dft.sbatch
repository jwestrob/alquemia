#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J GGR-static-context
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
ENVREP_MANIFEST=$1
case "$ENVREP_MANIFEST" in /*) ;; *) exit 2;; esac
ENVREP_OUT=$(dirname "$ENVREP_MANIFEST")
ENVREP_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
"$ENVREP_PY" "$ENVREP_OUT/implementation/affordable_workflow.py" execute --manifest "$ENVREP_MANIFEST" --output "$ENVREP_OUT/execution_${SLURM_JOB_ID}.json"
"$ENVREP_PY" "$ENVREP_OUT/implementation/environment_replica_context.py" collect --manifest "$ENVREP_MANIFEST" --output "$ENVREP_OUT/collection_${SLURM_JOB_ID}.json"
