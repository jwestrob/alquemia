#!/bin/bash
#SBATCH -p gpu
#SBATCH -J PQQ-context-H
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
SECOND_MANIFEST=$1
case "$SECOND_MANIFEST" in /*) ;; *) exit 2;; esac
SECOND_OUT=$(dirname "$SECOND_MANIFEST")
SECOND_PY=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1
"$SECOND_PY" "$SECOND_OUT/implementation/second_shell_hydrogen.py" execute --manifest "$SECOND_MANIFEST"
