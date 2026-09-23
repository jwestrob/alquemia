#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-lanm-series
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
lanm_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
exec "$lanm_root/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" \
  "$lanm_root/workspaces/lanm_series_followup_20260923/prepared_v2/implementation/lanm_series_followup.py" mace "$@"
