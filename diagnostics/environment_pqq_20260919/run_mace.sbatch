#!/bin/bash
#SBATCH -p gpu
#SBATCH -J PQQ-compact-context
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
ENVPQQ_MANIFEST=$1
case "$ENVPQQ_MANIFEST" in /*) ;; *) exit 2;; esac
ENVPQQ_OUT=$(dirname "$ENVPQQ_MANIFEST")
ENVPQQ_PY=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$ENVPQQ_PY" "$ENVPQQ_OUT/implementation/environment_pqq_context.py" execute --manifest "$ENVPQQ_MANIFEST"
"$ENVPQQ_PY" "$ENVPQQ_OUT/implementation/environment_pqq_context.py" collect --collection "$ENVPQQ_OUT/collection_${SLURM_JOB_ID}.json" --output "$ENVPQQ_OUT/result_${SLURM_JOB_ID}.json"
