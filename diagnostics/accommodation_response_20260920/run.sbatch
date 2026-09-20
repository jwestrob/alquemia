#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J accommodation-grad
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
ACCOMMODATION_MANIFEST=${1:?absolute manifest path required}
case "$ACCOMMODATION_MANIFEST" in /*) ;; *) exit 2;; esac
ACCOMMODATION_DIR=$(dirname "$ACCOMMODATION_MANIFEST")
ACCOMMODATION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
if "$ACCOMMODATION_PY" "$ACCOMMODATION_DIR/implementation/accommodation_response.py" execute --manifest "$ACCOMMODATION_MANIFEST"; then
    ACCOMMODATION_STATUS=0
else
    ACCOMMODATION_STATUS=$?
fi
"$ACCOMMODATION_PY" "$ACCOMMODATION_DIR/implementation/accommodation_response.py" collect \
    --manifest "$ACCOMMODATION_MANIFEST" --output "$ACCOMMODATION_DIR/collection_${SLURM_JOB_ID}.json"
exit "$ACCOMMODATION_STATUS"
