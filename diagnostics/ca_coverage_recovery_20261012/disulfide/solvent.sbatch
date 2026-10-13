#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=nikasha-MMOL2048-solvent
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --no-requeue
#SBATCH --output=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ca_coverage_recovery_20261012/disulfide/solvent_%j.out
#SBATCH --error=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ca_coverage_recovery_20261012/disulfide/solvent_%j.err
set -euo pipefail
test "$#" -eq 1
test -n "${SLURM_JOB_ID:-}"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python3 -u /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/ca_coverage_recovery_20261012/disulfide/SOLVENT.py --config "$1"
