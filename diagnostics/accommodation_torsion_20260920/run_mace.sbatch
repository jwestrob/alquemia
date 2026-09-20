#!/bin/bash
#SBATCH -p gpu
#SBATCH -J accommodation-torsion-MACE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
TORSION_MANIFEST=${1:?absolute manifest required}
case "$TORSION_MANIFEST" in /*) ;; *) exit 2;; esac
TORSION_ROOT=$(dirname "$TORSION_MANIFEST")
TORSION_PY=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$TORSION_PY" "$TORSION_ROOT/implementation/accommodation_torsion_profiles.py" mace --manifest "$TORSION_MANIFEST" --output "$TORSION_ROOT/mace_${SLURM_JOB_ID}"
