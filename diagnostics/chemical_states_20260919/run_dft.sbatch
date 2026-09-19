#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J chemical-proton
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: run_dft.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
CHEMICAL_MANIFEST=$1
case "$CHEMICAL_MANIFEST" in /*) ;; *) echo 'Manifest must be absolute' >&2; exit 2;; esac
CHEMICAL_OUTPUT=$(dirname "$CHEMICAL_MANIFEST")
CHEMICAL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
"$CHEMICAL_PY" scripts/affordable_workflow.py execute --manifest "$CHEMICAL_MANIFEST" --output "$CHEMICAL_OUTPUT/execution_${SLURM_JOB_ID}.json"
PYTHONPATH="$CHEMICAL_OUTPUT/implementation:$CHEMICAL_OUTPUT/implementation_${SLURM_JOB_ID}" "$CHEMICAL_PY" "$CHEMICAL_OUTPUT/implementation/chemical_states_proton.py" collect --manifest "$CHEMICAL_MANIFEST" --output "$CHEMICAL_OUTPUT/collection_${SLURM_JOB_ID}.json"
