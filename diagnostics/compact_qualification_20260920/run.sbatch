#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J compact-SCF-check
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
QUAL_MANIFEST=${1:?absolute manifest path required}
case "$QUAL_MANIFEST" in /*) ;; *) exit 2;; esac
QUAL_DIR=$(dirname "$QUAL_MANIFEST")
QUAL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
if "$QUAL_PY" "$QUAL_DIR/implementation/compact_solvation_qualification.py" execute --manifest "$QUAL_MANIFEST"; then
  QUAL_EXEC_STATUS=0
else
  QUAL_EXEC_STATUS=$?
fi
"$QUAL_PY" "$QUAL_DIR/implementation/compact_solvation_qualification.py" collect --manifest "$QUAL_MANIFEST" --output "$QUAL_DIR/collection_${SLURM_JOB_ID}.json"
exit "$QUAL_EXEC_STATUS"
