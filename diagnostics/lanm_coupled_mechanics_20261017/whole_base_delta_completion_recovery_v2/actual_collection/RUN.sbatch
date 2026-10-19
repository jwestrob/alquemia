#!/bin/bash
#SBATCH --account=banfield
#SBATCH --qos=standard
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH OMP_PLACES OMP_PROC_BIND
exec /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python "$@"
