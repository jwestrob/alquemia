#!/bin/bash
#SBATCH -p gpu
#SBATCH -J strict-donor-response
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
DONOR_MANIFEST=$1
case "$DONOR_MANIFEST" in /*) ;; *) exit 2;; esac
DONOR_ROOT=$(dirname "$DONOR_MANIFEST")
DONOR_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$DONOR_PY" "$DONOR_ROOT/implementation/strict_donor_response.py" execute --manifest "$DONOR_MANIFEST"
DONOR_STATUS=$?
"$DONOR_PY" "$DONOR_ROOT/implementation/strict_donor_response.py" collect --manifest "$DONOR_MANIFEST" --output "$DONOR_ROOT/COLLECTION.json"
DONOR_COLLECT=$?
if [ "$DONOR_STATUS" -ne 0 ]; then exit "$DONOR_STATUS"; fi
exit "$DONOR_COLLECT"
