#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-baseline
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/baseline_benchmark_20260915/run_v1/slurm_%j.out
#SBATCH -e /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/baseline_benchmark_20260915/run_v1/slurm_%j.err
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
MANIFEST=workspaces/baseline_benchmark_20260915/run_v1/manifest.json
"$PYTHON" scripts/affordable_workflow.py execute --manifest "$MANIFEST" --output "workspaces/baseline_benchmark_20260915/run_v1/execution_${SLURM_JOB_ID}.json"
"$PYTHON" scripts/affordable_benchmark.py collect --manifest "$MANIFEST" --output "workspaces/baseline_benchmark_20260915/run_v1/collection_${SLURM_JOB_ID}.json"
