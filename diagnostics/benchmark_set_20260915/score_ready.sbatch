#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-benchmark-set
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/benchmark_set_20260915/ready_tasks_v4/slurm_%j.out
#SBATCH -e /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/benchmark_set_20260915/ready_tasks_v4/slurm_%j.err
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
BENCHMARK_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
BENCHMARK_MANIFEST=workspaces/benchmark_set_20260915/ready_tasks_v4/manifest.json
"$BENCHMARK_PY" scripts/affordable_workflow.py execute --manifest "$BENCHMARK_MANIFEST" --output "workspaces/benchmark_set_20260915/ready_tasks_v4/execution_${SLURM_JOB_ID}.json"
"$BENCHMARK_PY" scripts/affordable_benchmark.py collect --manifest "$BENCHMARK_MANIFEST" --output "workspaces/benchmark_set_20260915/ready_tasks_v4/collection_${SLURM_JOB_ID}.json"
