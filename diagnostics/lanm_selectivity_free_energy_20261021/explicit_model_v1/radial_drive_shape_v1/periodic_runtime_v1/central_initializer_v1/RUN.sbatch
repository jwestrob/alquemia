#!/bin/bash
#SBATCH --partition=standard
#SBATCH --qos=standard_expand
#SBATCH --account=banfield
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
#SBATCH --job-name=periodic_v2_init
set -euo pipefail
unset LD_PRELOAD LD_LIBRARY_PATH
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
exec /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/_envs/lanm_cholmod_v1/bin/python "$1" --config "$2" --release "$3" --out "$4"
