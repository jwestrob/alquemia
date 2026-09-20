#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nonlinear-donors
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
NONLINEAR_MANIFEST=$1
case "$NONLINEAR_MANIFEST" in /*) ;; *) exit 2;; esac
NONLINEAR_RUN=$(dirname "$NONLINEAR_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python "$NONLINEAR_RUN/implementation/accommodation_nonlinear.py" execute --manifest "$NONLINEAR_MANIFEST"
