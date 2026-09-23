#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-union-mace
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
UNION_MANIFEST=$1
case "$UNION_MANIFEST" in /*) ;; *) exit 2;; esac
UNION_ROOT=$(dirname "$(dirname "$UNION_MANIFEST")")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNION_GPU_PY=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python
"$UNION_GPU_PY" "$UNION_ROOT/implementation/accommodation_folds.py" run-mace --manifest "$UNION_MANIFEST"
