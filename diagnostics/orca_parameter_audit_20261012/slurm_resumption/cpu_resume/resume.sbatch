#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J nikasha-ca-resume
#SBATCH --ntasks=14
#SBATCH --nodes=1
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/orca_parameter_audit_20261012/slurm_resumption_20261012/cpu_resume_%j.out

set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -u diagnostics/orca_parameter_audit_20261012/slurm_resumption/cpu_resume/LAUNCH.py --manifest workspaces/orca_parameter_audit_20261012/cpu_proposal_batch_preflight_v1/manifest.json --drain workspaces/orca_parameter_audit_20261012/slurm_resumption_20261012/DRAIN.json --receipts workspaces/orca_parameter_audit_20261012/slurm_resumption_20261012/cpu_resume
