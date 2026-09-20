#!/bin/bash
#SBATCH -p gpu
#SBATCH -J compact-score-timing
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
COMPACT_SCANNER_MANIFEST=$1
case "$COMPACT_SCANNER_MANIFEST" in /*) ;; *) exit 2;; esac
COMPACT_SCANNER_DIR=$(dirname "$COMPACT_SCANNER_MANIFEST")
COMPACT_SCANNER_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$COMPACT_SCANNER_PY" "$COMPACT_SCANNER_DIR/implementation/compact_solvation_scanner.py" execute --manifest "$COMPACT_SCANNER_MANIFEST"
