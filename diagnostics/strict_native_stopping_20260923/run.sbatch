#!/bin/bash
#SBATCH -p gpu
#SBATCH -J strict-native-stopping
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
STRICT_MANIFEST=$1
case "$STRICT_MANIFEST" in /*) ;; *) exit 2;; esac
STRICT_ROOT=$(dirname "$STRICT_MANIFEST")
STRICT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
STRICT_SCRIPT="$STRICT_ROOT/implementation/strict_native_stopping.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$STRICT_PY" "$STRICT_SCRIPT" execute --manifest "$STRICT_MANIFEST"
STRICT_EXEC=$?
"$STRICT_PY" "$STRICT_SCRIPT" collect --manifest "$STRICT_MANIFEST" --output "$STRICT_ROOT/COLLECTION.json"
STRICT_COLLECT=$?
if [ "$STRICT_EXEC" -ne 0 ]; then exit "$STRICT_EXEC"; fi
exit "$STRICT_COLLECT"
