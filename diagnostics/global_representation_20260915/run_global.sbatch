#!/bin/bash
#SBATCH -p standard,memory,high-memory
#SBATCH -J alquemia-global-xtb
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/global_representation_20260915/slurm_%j.out
#SBATCH -e /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/global_representation_20260915/slurm_%j.err
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
GLOBAL_XTB_RANKS=$((SLURM_CPUS_ON_NODE / 2))
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/global_representation_20260915/implementation/run_orca_task_manifest.py /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/global_representation_20260915/prepared_v1/global_manifest.json --orca /groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca --workers 2 --nprocs "$GLOBAL_XTB_RANKS"
