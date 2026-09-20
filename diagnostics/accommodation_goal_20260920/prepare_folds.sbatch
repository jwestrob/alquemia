#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J PQQ-fold-prepare
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
set -euo pipefail
FOLD_RUN=${1:?absolute run directory}
FOLD_ROOT=${2:?absolute repository}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
FOLD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$FOLD_PY" "$FOLD_RUN/implementation/accommodation_folds.py" prepare --sources "$FOLD_ROOT/diagnostics/accommodation_controls_20260920/PQQ_ALL250_SOURCES.json" --agreement "$FOLD_ROOT/diagnostics/accommodation_goal_20260920/FOLD_ROBUSTNESS_PLAN.md" --output "$FOLD_RUN/source_preparation" --workers 16
"$FOLD_PY" "$FOLD_RUN/implementation/accommodation_folds.py" prepare-mace --preparation "$FOLD_RUN/source_preparation/preparation.json" --agreement "$FOLD_ROOT/diagnostics/accommodation_goal_20260920/FOLD_ROBUSTNESS_PLAN.md" --output "$FOLD_RUN/mace"
