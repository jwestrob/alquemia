#!/bin/bash
#SBATCH -p gpu
#SBATCH -J triple-q0-mace
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
TRIPLE_STAGE=$1
case "$TRIPLE_STAGE" in /*) ;; *) exit 2;; esac
TRIPLE_GPU=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$TRIPLE_GPU" "$TRIPLE_STAGE/implementation/accommodation_folds.py" run-mace --manifest "$TRIPLE_STAGE/mace/manifest.json"
