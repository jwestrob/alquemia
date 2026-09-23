#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-ftol-pilot
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
PRECISION_MANIFEST=$1
case "$PRECISION_MANIFEST" in /*) ;; *) exit 2;; esac
PRECISION_RUN=$(dirname "$PRECISION_MANIFEST")
PRECISION_ROOT=$(dirname "$PRECISION_RUN")
PRECISION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
PRECISION_SCRIPT="$PRECISION_RUN/implementation/slsqp_precision.py"
PRECISION_AGREEMENT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/slsqp_precision_20260923/PLAN.md
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$PRECISION_PY" "$PRECISION_SCRIPT" execute --manifest "$PRECISION_MANIFEST"
"$PRECISION_PY" "$PRECISION_SCRIPT" collect --manifest "$PRECISION_MANIFEST" --output "$PRECISION_RUN/collection_${SLURM_JOB_ID}.json"
"$PRECISION_PY" "$PRECISION_SCRIPT" prepare_pool --proposals "$PRECISION_RUN/collection_${SLURM_JOB_ID}.json" --agreement "$PRECISION_AGREEMENT" --output "$PRECISION_ROOT/pool_v1"
"$PRECISION_PY" "$PRECISION_ROOT/pool_v1/implementation/slsqp_precision.py" execute_mace --manifest "$PRECISION_ROOT/pool_v1/manifest.json"
