#!/bin/bash
#SBATCH --job-name=lady-boundary-one
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --qos=standard_expand
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=13
#SBATCH --no-requeue
set -euo pipefail
unset LD_PRELOAD LD_LIBRARY_PATH
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
test "$#" -eq 2
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
exec workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/conic_saved_comparison_v1/environment_v1/venv/bin/python diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/radial_drive_shape_v1/boundary_fraction_witness_v1/RUN_v2.py run --config "$1" --out "$2"
