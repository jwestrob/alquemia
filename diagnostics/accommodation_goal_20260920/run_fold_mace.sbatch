#!/bin/bash
#SBATCH -p gpu
#SBATCH -J PQQ-fold-MACE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
FOLD_RUN=${1:?absolute run directory}
FOLD_ROOT=${2:?absolute repository}
FOLD_MACE=${3:?absolute MACE directory}
FOLD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FOLD_GPU_PY="$FOLD_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$FOLD_GPU_PY" "$FOLD_MACE/implementation/accommodation_folds.py" run-mace --manifest "$FOLD_MACE/manifest.json"
"$FOLD_PY" "$FOLD_MACE/implementation/accommodation_folds.py" inventory --manifest "$FOLD_MACE/manifest.json" --output "$FOLD_RUN/INVENTORY.json"
"$FOLD_PY" "$FOLD_ROOT/scripts/compact_solvation.py" prepare --inventory "$FOLD_RUN/INVENTORY.json" --agreement "$FOLD_ROOT/diagnostics/accommodation_goal_20260920/FOLD_ROBUSTNESS_PLAN.md" --output "$FOLD_RUN/solvent" --reuse-manifest "$FOLD_ROOT/workspaces/compact_solvation_20260920/full_v1/manifest.json"
