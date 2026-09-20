#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J accommodation-torsion
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
TORSION_MANIFEST=${1:?absolute stage manifest required}
case "$TORSION_MANIFEST" in /*) ;; *) exit 2;; esac
TORSION_STAGE=$(dirname "$TORSION_MANIFEST")
TORSION_ROOT=$(dirname "$TORSION_STAGE")
TORSION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
if "$TORSION_PY" "$TORSION_ROOT/implementation/accommodation_torsion_profiles.py" execute --manifest "$TORSION_MANIFEST"; then
 TORSION_STATUS=0
else
 TORSION_STATUS=$?
fi
"$TORSION_PY" "$TORSION_ROOT/implementation/accommodation_torsion_profiles.py" collect_native --manifest "$TORSION_MANIFEST" --output "$TORSION_STAGE/collection_${SLURM_JOB_ID}.json"
exit "$TORSION_STATUS"
