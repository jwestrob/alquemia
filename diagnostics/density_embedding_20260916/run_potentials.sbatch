#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-density-coupling
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
set -euo pipefail
if [ "$#" -ne 2 ]; then
    echo 'Usage: run_potentials.sbatch ABSOLUTE_MANIFEST ABSOLUTE_OUTPUT_DIR' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
DENSITY_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$DENSITY_PY" scripts/density_embedding.py execute-potentials --manifest "$1" --output "$2/execution_${SLURM_JOB_ID}.json"
"$DENSITY_PY" scripts/density_embedding.py collect-potentials --manifest "$1" --output "$2/potentials_${SLURM_JOB_ID}.json"
