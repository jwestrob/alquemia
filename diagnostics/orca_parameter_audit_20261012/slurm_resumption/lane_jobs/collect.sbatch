#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J nikasha-ca-collect
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --nodelist=node-64-768g-1
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/orca_parameter_audit_20261012/slurm_resumption_20261012/collect_%j.out
set -euo pipefail
test "$#" -eq 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -u diagnostics/orca_parameter_audit_20261012/slurm_resumption/lane_jobs/LANES.py collect --plan "$1"
