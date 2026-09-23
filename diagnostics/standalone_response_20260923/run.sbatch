#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-solvent-response
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
RESPONSE_MANIFEST=$1
case "$RESPONSE_MANIFEST" in /*) ;; *) exit 2;; esac
RESPONSE_RUN=$(dirname "$RESPONSE_MANIFEST")
RESPONSE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
RESPONSE_SCRIPT="$RESPONSE_RUN/implementation/standalone_response_proposals.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$RESPONSE_PY" "$RESPONSE_SCRIPT" execute --manifest "$RESPONSE_MANIFEST"
RESPONSE_STATUS=$?
if [ -f "$RESPONSE_RUN/collection.json" ]; then
 "$RESPONSE_PY" "$RESPONSE_SCRIPT" report --collection "$RESPONSE_RUN/collection.json" --output "$RESPONSE_RUN/COMPARISON.json"
fi
exit "$RESPONSE_STATUS"
