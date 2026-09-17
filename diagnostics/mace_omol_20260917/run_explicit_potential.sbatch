#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-explicit-potential
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: run_explicit_potential.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FIELD_DIR=$(dirname "$1")
"$FIELD_PY" "$FIELD_DIR/implementation/mace_explicit_field_short.py" execute-potential --manifest "$1" --output "$FIELD_DIR/potential_execution_${SLURM_JOB_ID}.json"
