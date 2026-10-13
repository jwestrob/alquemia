#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=Nikasha-MMOL2048-recovery
#SBATCH --nodes=1
#SBATCH --ntasks=7
#SBATCH --no-requeue
#SBATCH --output=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ca_coverage_recovery_20261012/disulfide/native_v2/slurm_%j.out
#SBATCH --error=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ca_coverage_recovery_20261012/disulfide/native_v2/slurm_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Scientific computation requires Slurm}"
if [[ $# -ne 1 ]]; then echo 'Usage: run_native.sbatch ABSOLUTE_MANIFEST_JSON' >&2; exit 2; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/ca_coverage_recovery_20261012/disulfide/NATIVE.py execute --manifest "$1"
