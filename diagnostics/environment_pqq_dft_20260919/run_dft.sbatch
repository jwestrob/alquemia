#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J PQQ-context-charge
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
ENVPQQ_DFT_MANIFEST=$1
case "$ENVPQQ_DFT_MANIFEST" in /*) ;; *) exit 2;; esac
ENVPQQ_DFT_OUT=$(dirname "$ENVPQQ_DFT_MANIFEST")
ENVPQQ_DFT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
"$ENVPQQ_DFT_PY" "$ENVPQQ_DFT_OUT/implementation/affordable_workflow.py" execute --manifest "$ENVPQQ_DFT_MANIFEST" --output "$ENVPQQ_DFT_OUT/execution_${SLURM_JOB_ID}.json"
"$ENVPQQ_DFT_PY" "$ENVPQQ_DFT_OUT/implementation/environment_pqq_dft.py" collect --manifest "$ENVPQQ_DFT_MANIFEST" --output "$ENVPQQ_DFT_OUT/collection_${SLURM_JOB_ID}.json"
