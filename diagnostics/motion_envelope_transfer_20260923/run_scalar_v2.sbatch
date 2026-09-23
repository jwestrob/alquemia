#!/bin/bash
#SBATCH -p gpu,standard-shared
#SBATCH -J envelope100-strict
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=65536M
set -euo pipefail
ENVELOPE_REQUESTED_IMPL=$1
ENVELOPE_MANIFEST=$2
ENVELOPE_IMPL="${ENVELOPE_MANIFEST%/*}/implementation"
case "$ENVELOPE_IMPL:$ENVELOPE_MANIFEST" in /*:/*) ;; *) exit 2;; esac
ENVELOPE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
set +e
"$ENVELOPE_CPU" "$ENVELOPE_IMPL/motion_envelope_transfer.py" execute_scalar --manifest "$ENVELOPE_MANIFEST"
ENVELOPE_RC=$?
set -e
"$ENVELOPE_CPU" "$ENVELOPE_IMPL/motion_envelope_transfer.py" collect_scalar --manifest "$ENVELOPE_MANIFEST" --output "${ENVELOPE_MANIFEST%/*}/COLLECTION.json"
exit "$ENVELOPE_RC"
