#!/bin/bash
#SBATCH -p gpu,standard-shared
#SBATCH -J envelope-strict
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=65536M
set -euo pipefail
ENVELOPE_STAGE=$1
case "$ENVELOPE_STAGE" in /*) ;; *) exit 2;; esac
ENVELOPE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
set +e
"$ENVELOPE_CPU" "$ENVELOPE_STAGE/implementation/motion_envelope_scalar.py" execute --manifest "$ENVELOPE_STAGE/manifest.json"
ENVELOPE_RC=$?
set -e
"$ENVELOPE_CPU" "$ENVELOPE_STAGE/implementation/motion_envelope_scalar.py" collect --manifest "$ENVELOPE_STAGE/manifest.json" --output "$ENVELOPE_STAGE/COLLECTION.json"
exit "$ENVELOPE_RC"
