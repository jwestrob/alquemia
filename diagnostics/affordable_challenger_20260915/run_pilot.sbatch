#!/bin/bash
#SBATCH -p standard,memory,high-memory
#SBATCH -J alquemia-mbis-pilot
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/affordable_challenger_20260915/pilot/slurm_%j.out
#SBATCH -e /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/affordable_challenger_20260915/pilot/slurm_%j.err
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
MANIFEST=workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json
"$PYTHON" scripts/affordable_workflow.py execute --manifest "$MANIFEST" --output "workspaces/affordable_challenger_20260915/pilot/execution_${SLURM_JOB_ID}.json"
"$PYTHON" scripts/affordable_workflow.py collect --manifest "$MANIFEST" --output "workspaces/affordable_challenger_20260915/pilot/collection_${SLURM_JOB_ID}.json"
