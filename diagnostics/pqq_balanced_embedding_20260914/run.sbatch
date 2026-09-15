#!/bin/bash
#SBATCH -p high-memory
#SBATCH -J pqq-balanced-PC
#SBATCH --exclusive
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_balanced_embedding_20260914/slurm_%j.out
#SBATCH -e /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_balanced_embedding_20260914/slurm_%j.err

set -euo pipefail
ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_balanced_embedding_20260914
PROJECT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
ORCA_PATH=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export ORCA_PATH
export PATH=$ORCA_PATH:$PATH
export LD_LIBRARY_PATH=$ORCA_PATH:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$ROOT"

"$PYTHON" prepare_embedding.py verify
WORKERS=12
NPROCS=$((SLURM_CPUS_ON_NODE / WORKERS))
if (( NPROCS < 1 )); then
    echo "allocation too small for the 12-task panel" >&2
    exit 2
fi
"$PYTHON" "$PROJECT/scripts/run_orca_task_manifest.py" embedded_tasks.json \
    --orca "$ORCA_PATH/orca" --workers "$WORKERS" --nprocs "$NPROCS" \
    --expected-runner-sha256 b9278e74e5317e859e8cb05541e1ec57cb335b26dcdf65a83b89675407343b9a
"$PYTHON" prepare_embedding.py verify
"$PYTHON" summarize.py
