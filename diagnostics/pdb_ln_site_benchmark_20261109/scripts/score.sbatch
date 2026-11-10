#!/bin/bash
#SBATCH -J lnbench_score
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard
#SBATCH -c 12
#SBATCH -N 1
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pdb_ln_site_benchmark_20261109/logs/score_%j.out

set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
B=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
D=$B/diagnostics/pdb_ln_site_benchmark_20261109
W=$B/workspaces/pdb_ln_site_benchmark_20261109
cd $W
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -I $D/scripts/score.py $W $D/scripts $B/workspaces/laca_imprint_audit_20261109/cif /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner $SLURM_CPUS_ON_NODE
