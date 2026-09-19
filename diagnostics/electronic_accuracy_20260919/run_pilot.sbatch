#!/bin/bash
#SBATCH --job-name=electronic_cc_v1
#SBATCH --partition=standard,memory
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
MANIFEST=${1:?Provide absolute manifest path}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
cd "$ROOT"
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/affordable_workflow.py execute --manifest "$MANIFEST" --output "$(dirname "$MANIFEST")/execution_${SLURM_JOB_ID}.json"
