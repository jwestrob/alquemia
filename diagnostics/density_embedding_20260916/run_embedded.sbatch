#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-density-response
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [ "$#" -ne 2 ]; then
    echo 'Usage: run_embedded.sbatch ABSOLUTE_MANIFEST ABSOLUTE_OUTPUT_DIR' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
DENSITY_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$DENSITY_PY" scripts/density_embedding.py execute-embedded --manifest "$1" --output "$2/execution_${SLURM_JOB_ID}.json"
