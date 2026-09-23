#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-six-strict
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
SIX_MANIFEST=$1
case "$SIX_MANIFEST" in /*) ;; *) exit 2;; esac
SIX_RUN=$(dirname "$SIX_MANIFEST")
SIX_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
SIX_SCRIPT="$SIX_RUN/implementation/six_angle_accommodation.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$SIX_PY" "$SIX_SCRIPT" execute_low --manifest "$SIX_MANIFEST"
SIX_STATUS=$?
"$SIX_PY" "$SIX_SCRIPT" collect --manifest "$SIX_MANIFEST" --output "$SIX_RUN/COLLECTION.json"
SIX_COLLECT=$?
if [ "$SIX_STATUS" -ne 0 ]; then exit "$SIX_STATUS"; fi
exit "$SIX_COLLECT"
