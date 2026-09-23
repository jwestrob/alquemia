#!/bin/bash
#SBATCH -p gpu
#SBATCH -J strict-native225
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
set -uo pipefail
STRICT_MANIFEST=$1
case "$STRICT_MANIFEST" in /*) ;; *) exit 2;; esac
STRICT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
STRICT_SCRIPT=$($STRICT_PY -c 'import json,sys; print(json.load(open(sys.argv[1]))["implementation"]["strict_native_transfer.py"]["path"])' "$STRICT_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
STRICT_STATUS=0
"$STRICT_PY" "$STRICT_SCRIPT" execute --manifest "$STRICT_MANIFEST" || STRICT_STATUS=$?
"$STRICT_PY" "$STRICT_SCRIPT" collect --manifest "$STRICT_MANIFEST" --output "${STRICT_MANIFEST%/*}/COLLECTION.json" || { if [ "$STRICT_STATUS" -eq 0 ]; then STRICT_STATUS=1; fi; }
exit "$STRICT_STATUS"
