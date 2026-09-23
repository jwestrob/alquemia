#!/bin/bash
#SBATCH -p gpu
#SBATCH -J strict-native-pool32
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
STRICT_ROOT=$1
case "$STRICT_ROOT" in /*) ;; *) exit 2;; esac
STRICT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
STRICT_STATUS=0
for STRICT_BRANCH in fresh cold_seed; do
  STRICT_MANIFEST="$STRICT_ROOT/$STRICT_BRANCH/manifest.json"
  STRICT_SCRIPT="$STRICT_ROOT/$STRICT_BRANCH/implementation/strict_native_pool.py"
  "$STRICT_PY" "$STRICT_SCRIPT" execute --manifest "$STRICT_MANIFEST" || STRICT_STATUS=1
  "$STRICT_PY" "$STRICT_SCRIPT" collect --manifest "$STRICT_MANIFEST" --output "$STRICT_ROOT/$STRICT_BRANCH/COLLECTION.json" || STRICT_STATUS=1
done
exit "$STRICT_STATUS"
