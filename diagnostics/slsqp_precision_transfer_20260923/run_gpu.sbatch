#!/bin/bash
#SBATCH -p gpu
#SBATCH -J precision225
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
TRANSFER_MANIFEST=$1
case "$TRANSFER_MANIFEST" in /*) ;; *) exit 2;; esac
TRANSFER_RUN=$(dirname "$TRANSFER_MANIFEST")
TRANSFER_ROOT=$(dirname "$TRANSFER_RUN")
TRANSFER_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
TRANSFER_SCRIPT="$TRANSFER_RUN/implementation/slsqp_precision_transfer.py"
TRANSFER_AGREEMENT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/slsqp_precision_transfer_20260923/PLAN.md
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$TRANSFER_PY" "$TRANSFER_SCRIPT" execute --manifest "$TRANSFER_MANIFEST"
"$TRANSFER_PY" "$TRANSFER_SCRIPT" collect --manifest "$TRANSFER_MANIFEST" --output "$TRANSFER_RUN/collection_final.json"
"$TRANSFER_PY" "$TRANSFER_SCRIPT" prepare_pool --proposals "$TRANSFER_RUN/collection_final.json" --agreement "$TRANSFER_AGREEMENT" --output "$TRANSFER_ROOT/pool"
"$TRANSFER_PY" "$TRANSFER_ROOT/pool/implementation/slsqp_precision_transfer.py" execute_mace --manifest "$TRANSFER_ROOT/pool/manifest.json"
