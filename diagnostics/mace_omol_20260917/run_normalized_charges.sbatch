#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-normalized-charges
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: run_normalized_charges.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
CHARGE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
CHARGE_DIR=$(dirname "$1")
"$CHARGE_PY" "$CHARGE_DIR/implementation/mace_omol_charges.py" execute --manifest "$1" --output "$CHARGE_DIR/execution_${SLURM_JOB_ID}.json"
