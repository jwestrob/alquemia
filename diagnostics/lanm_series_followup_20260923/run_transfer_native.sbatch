#!/bin/bash
#SBATCH -p gpu
#SBATCH -J lanm-Dy-transfer-gfn
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
lanm_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 "$lanm_root/workspaces/lanm_series_followup_20260923/dy_transfer_v1/implementation/lanm_structural_transfer.py" execute-native "$@"
