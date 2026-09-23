#!/bin/bash
#SBATCH -p gpu
#SBATCH -J precision-self-continuation
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
CONT_MANIFEST=$1
case "$CONT_MANIFEST" in /*) ;; *) exit 2;; esac
CONT_STAGE1=$(dirname "$CONT_MANIFEST")
CONT_ROOT=$(dirname "$CONT_STAGE1")
CONT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
CONT_SCRIPT="$CONT_STAGE1/implementation/precision_geometry_continuation.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$CONT_PY" "$CONT_SCRIPT" execute --manifest "$CONT_MANIFEST"
CONT_EXEC1=$?
"$CONT_PY" "$CONT_SCRIPT" collect --manifest "$CONT_MANIFEST" --output "$CONT_STAGE1/collection.json"
CONT_COLLECT1=$?
if [ "$CONT_COLLECT1" -ne 0 ]; then exit "$CONT_COLLECT1"; fi
"$CONT_PY" "$CONT_SCRIPT" prepare --sources "$2" --agreement "$3" --previous "$CONT_STAGE1/collection.json" --output "$CONT_ROOT/stage2"
CONT_PREP2=$?
if [ "$CONT_PREP2" -ne 0 ]; then exit "$CONT_PREP2"; fi
CONT_SCRIPT2="$CONT_ROOT/stage2/implementation/precision_geometry_continuation.py"
"$CONT_PY" "$CONT_SCRIPT2" execute --manifest "$CONT_ROOT/stage2/manifest.json"
CONT_EXEC2=$?
"$CONT_PY" "$CONT_SCRIPT2" collect --manifest "$CONT_ROOT/stage2/manifest.json" --output "$CONT_ROOT/stage2/collection.json"
CONT_COLLECT2=$?
if [ "$CONT_COLLECT2" -ne 0 ]; then exit "$CONT_COLLECT2"; fi
"$CONT_PY" "$CONT_SCRIPT2" compare --stage1 "$CONT_STAGE1/collection.json" --stage2 "$CONT_ROOT/stage2/collection.json" --output "$CONT_ROOT/COMPARISON.json"
CONT_COMPARE=$?
if [ "$CONT_EXEC1" -ne 0 ]; then exit "$CONT_EXEC1"; fi
if [ "$CONT_EXEC2" -ne 0 ]; then exit "$CONT_EXEC2"; fi
exit "$CONT_COMPARE"
