#!/bin/bash
#SBATCH -p gpu,standard-shared
#SBATCH -J envelope100-collect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192M
set -euo pipefail
ENVELOPE_POOL=$1
case "$ENVELOPE_POOL" in /*) ;; *) exit 2;; esac
ENVELOPE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$ENVELOPE_CPU" "$ENVELOPE_POOL/implementation/motion_envelope_transfer.py" collect_pool --manifest "$ENVELOPE_POOL/manifest.json" --output "$ENVELOPE_POOL/COLLECTION_FINAL.json"
