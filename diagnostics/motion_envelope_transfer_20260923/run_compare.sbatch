#!/bin/bash
#SBATCH -p gpu,standard-shared
#SBATCH -J envelope100-compare
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192M
set -euo pipefail
ENVELOPE_ROOT=$1
case "$ENVELOPE_ROOT" in /*) ;; *) exit 2;; esac
ENVELOPE_REPO=$2
ENVELOPE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$ENVELOPE_CPU" "$ENVELOPE_ROOT/analysis_v1/implementation/motion_envelope_transfer_compare.py" \
 --inventory "$ENVELOPE_ROOT/run_v1/INVENTORY.json" \
 --collections "$ENVELOPE_ROOT/run_v2/shard_0/pool/COLLECTION_FINAL.json" "$ENVELOPE_ROOT/run_v2/shard_1/pool/COLLECTION_FINAL.json" \
 --strict "$ENVELOPE_REPO/workspaces/strict_native_transfer_20260923/run_v1/COMPARISON.json" \
 --threefold "$ENVELOPE_REPO/workspaces/union_triple_transfer_20260923/COMPARISON_v1.json" \
 --output "$ENVELOPE_ROOT/COMPARISON_v1.json"
