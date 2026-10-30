#!/bin/bash
#SBATCH -A banfield
#SBATCH -p standard
#SBATCH --qos=standard_expand
#SBATCH -c 24
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=$SLURM_CPUS_ON_NODE PYTHONUNBUFFERED=1
export LADY_REPO_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
echo "node $(hostname) cpus $SLURM_CPUS_ON_NODE mem $SLURM_MEM_PER_NODE"
exec /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lady_almo_implementation_20261023/venv/bin/python "$@"
