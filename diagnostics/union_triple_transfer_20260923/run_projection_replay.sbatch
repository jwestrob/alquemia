#!/bin/bash
#SBATCH -p gpu
#SBATCH -J triple100-metadata
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192M
#SBATCH --nodelist=node-128-512g-8gpu-1
set -euo pipefail
TRIPLE_ROOT=$1
TRIPLE_DIAGNOSTIC=$2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$TRIPLE_ROOT/execution_v1:$TRIPLE_ROOT/searches/shard_0/proposals/implementation"
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python "$TRIPLE_DIAGNOSTIC" \
 "$TRIPLE_ROOT/SELECTOR_HOST_DIAGNOSTIC.json" \
 "$TRIPLE_ROOT/searches/shard_0/proposals/manifest.json" \
 "$TRIPLE_ROOT/searches/shard_1/proposals/manifest.json"
