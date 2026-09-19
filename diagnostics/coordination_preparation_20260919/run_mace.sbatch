#!/bin/bash
#SBATCH -p gpu
#SBATCH -J coordination-context
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
COORD_MANIFEST=$1
case "$COORD_MANIFEST" in /*) ;; *) exit 2;; esac
COORD_OUT=$(dirname "$COORD_MANIFEST")
COORD_PY=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$COORD_PY" "$COORD_OUT/implementation/coordination_preparation_context.py" execute --manifest "$COORD_MANIFEST"
