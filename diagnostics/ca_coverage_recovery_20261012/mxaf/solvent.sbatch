#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=nikasha-mxaf-solvent
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --no-requeue
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ca_coverage_recovery_20261012/mxaf_solvent_%j.out
set -euo pipefail
test "$#" -eq 4
test -n "${SLURM_JOB_ID:-}"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
python3 -u diagnostics/ca_coverage_recovery_20261012/mxaf/SOLVENT.py --scope "$1" --collection "$2" --primary "$3" --directory "$4"
