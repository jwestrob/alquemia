#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-pool-continuation
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
POOL_MANIFEST=$1
case "$POOL_MANIFEST" in /*) ;; *) exit 2;; esac
POOL_STAGE1=$(dirname "$POOL_MANIFEST")
POOL_ROOT=$(dirname "$POOL_STAGE1")
POOL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
POOL_SCRIPT="$POOL_STAGE1/implementation/native_pool_continuation.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$POOL_PY" "$POOL_SCRIPT" execute --manifest "$POOL_MANIFEST"
POOL_EXEC1=$?
"$POOL_PY" "$POOL_SCRIPT" collect --manifest "$POOL_MANIFEST" --output "$POOL_STAGE1/collection.json"
POOL_COLLECT1=$?
if [ "$POOL_COLLECT1" -ne 0 ]; then exit "$POOL_COLLECT1"; fi
"$POOL_PY" "$POOL_SCRIPT" prepare --inventory "$2" --agreement "$3" --previous "$POOL_STAGE1/collection.json" --output "$POOL_ROOT/stage2"
POOL_PREP2=$?
if [ "$POOL_PREP2" -ne 0 ]; then exit "$POOL_PREP2"; fi
POOL_SCRIPT2="$POOL_ROOT/stage2/implementation/native_pool_continuation.py"
"$POOL_PY" "$POOL_SCRIPT2" execute --manifest "$POOL_ROOT/stage2/manifest.json"
POOL_EXEC2=$?
"$POOL_PY" "$POOL_SCRIPT2" collect --manifest "$POOL_ROOT/stage2/manifest.json" --output "$POOL_ROOT/stage2/collection.json"
POOL_COLLECT2=$?
if [ "$POOL_COLLECT2" -ne 0 ]; then exit "$POOL_COLLECT2"; fi
"$POOL_PY" "$POOL_SCRIPT2" compare --stage1 "$POOL_STAGE1/collection.json" --stage2 "$POOL_ROOT/stage2/collection.json" --output "$POOL_ROOT/result.json"
POOL_COMPARE=$?
if [ "$POOL_EXEC1" -ne 0 ]; then exit "$POOL_EXEC1"; fi
if [ "$POOL_EXEC2" -ne 0 ]; then exit "$POOL_EXEC2"; fi
exit "$POOL_COMPARE"
