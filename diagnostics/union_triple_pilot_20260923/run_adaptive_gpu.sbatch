#!/bin/bash
#SBATCH -p gpu
#SBATCH -J triple-adaptive
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
TRIPLE_MANIFEST=$1
case "$TRIPLE_MANIFEST" in /*) ;; *) exit 2;; esac
TRIPLE_RUN=$(dirname "$TRIPLE_MANIFEST")
TRIPLE_ROOT=$(dirname "$TRIPLE_RUN")
TRIPLE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
TRIPLE_SCRIPT="$TRIPLE_RUN/implementation/union_triple_adaptive.py"
TRIPLE_PLAN=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/union_triple_pilot_20260923/PLAN.md
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" execute --manifest "$TRIPLE_MANIFEST"
TRIPLE_EXEC=$?
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" collect --manifest "$TRIPLE_MANIFEST" --output "$TRIPLE_RUN/collection.json" || exit $?
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" prepare_pool --proposals "$TRIPLE_RUN/collection.json" --agreement "$TRIPLE_PLAN" --output "$TRIPLE_ROOT/pool_v1" || exit $?
"$TRIPLE_CPU" "$TRIPLE_ROOT/pool_v1/implementation/union_triple_adaptive.py" execute_mace --manifest "$TRIPLE_ROOT/pool_v1/manifest.json"
TRIPLE_CROSS=$?
if [ "$TRIPLE_EXEC" -ne 0 ]; then exit "$TRIPLE_EXEC"; fi
exit "$TRIPLE_CROSS"
