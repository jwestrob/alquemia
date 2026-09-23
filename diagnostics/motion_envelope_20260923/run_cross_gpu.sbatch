#!/bin/bash
#SBATCH -p gpu
#SBATCH -J envelope-cross
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
ENVELOPE_STAGE=$1
case "$ENVELOPE_STAGE" in /*) ;; *) exit 2;; esac
ENVELOPE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$ENVELOPE_CPU" "$ENVELOPE_STAGE/implementation/motion_envelope_pool.py" execute_mace --manifest "$ENVELOPE_STAGE/manifest.json"
