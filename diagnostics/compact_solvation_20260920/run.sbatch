#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J compact-GFN2-solv
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
COMPACT_MANIFEST=${1:?absolute manifest path required}
case "$COMPACT_MANIFEST" in /*) ;; *) exit 2;; esac
COMPACT_DIR=$(dirname "$COMPACT_MANIFEST")
COMPACT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
"$COMPACT_PY" "$COMPACT_DIR/implementation/compact_solvation.py" execute --manifest "$COMPACT_MANIFEST"
"$COMPACT_PY" "$COMPACT_DIR/implementation/compact_solvation.py" collect --manifest "$COMPACT_MANIFEST" --output "$COMPACT_DIR/collection_${SLURM_JOB_ID}.json"
