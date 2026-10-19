#!/bin/bash
#SBATCH --job-name=lbt_base_delta_preflight
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH OMP_PLACES OMP_PROC_BIND
exec "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_coupled_mechanics_20261017/whole_base_delta_v1/BASE_DELTA.py" --stage "${1:-preflight}" --config "${2:-/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_coupled_mechanics_20261017/whole_base_delta_preflight_v1/CONFIG.json}" --output "${3:-/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_coupled_mechanics_20261017/whole_base_delta_preflight_v1/execution}"
