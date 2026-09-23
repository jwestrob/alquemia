#!/bin/bash
#SBATCH -p gpu
#SBATCH -J triple-q0-gfn
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=65536M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
TRIPLE_STAGE=$1
case "$TRIPLE_STAGE" in /*) ;; *) exit 2;; esac
TRIPLE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$TRIPLE_CPU" "$TRIPLE_STAGE/implementation/union_triple_origins.py" execute --manifest "$TRIPLE_STAGE/solvent/manifest.json"
