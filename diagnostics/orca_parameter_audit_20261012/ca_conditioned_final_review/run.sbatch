#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=Nikasha-Ca-final-review
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
#SBATCH --output=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/orca_parameter_audit_20261012/ca_conditioned_final_review_v1/slurm_%j.out
#SBATCH --error=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/orca_parameter_audit_20261012/ca_conditioned_final_review_v1/slurm_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Review requires Slurm}"
if [[ $# -ne 1 ]]; then echo 'Usage: run.sbatch ABSOLUTE_PLAN_JSON' >&2; exit 2; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python3 /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/orca_parameter_audit_20261012/ca_conditioned_final_review/REVIEW.py --plan "$1"
