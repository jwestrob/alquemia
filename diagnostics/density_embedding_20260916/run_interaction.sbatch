#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-asp-interaction
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: run_interaction.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
INTERACTION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
INTERACTION_DIR=$(dirname "$1")
"$INTERACTION_PY" scripts/affordable_workflow.py execute --manifest "$1" --output "$INTERACTION_DIR/execution_${SLURM_JOB_ID}.json"
