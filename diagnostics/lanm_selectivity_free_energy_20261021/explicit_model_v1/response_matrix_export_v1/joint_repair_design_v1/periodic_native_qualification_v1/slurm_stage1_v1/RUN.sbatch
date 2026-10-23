#!/bin/bash
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --qos=standard_expand
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
[[ "$#" -eq 3 ]]
exec env -u LD_PRELOAD -u LD_LIBRARY_PATH OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
  "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/_envs/lanm_cholmod_v1/bin/python" "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/periodic_native_qualification_v1/slurm_stage1_v1/RUN.py" --config "$1" --release "$2" --out "$3"
