#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-H-repair
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
COVERAGE_MANIFEST=$1
COVERAGE_PHASE=$2
case "$COVERAGE_MANIFEST" in /*) ;; *) exit 2;; esac
COVERAGE_RUN=$(dirname "$COVERAGE_MANIFEST")
COVERAGE_ROOT=$(dirname "$COVERAGE_RUN")
COVERAGE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
COVERAGE_SCRIPT="$COVERAGE_RUN/implementation/preparation_coverage_score.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
if [ "$COVERAGE_PHASE" = origins ]; then
  "$COVERAGE_PY" "$COVERAGE_SCRIPT" execute_origins --manifest "$COVERAGE_MANIFEST"
elif [ "$COVERAGE_PHASE" = proposals ]; then
  "$COVERAGE_PY" "$COVERAGE_SCRIPT" execute --manifest "$COVERAGE_MANIFEST"
  "$COVERAGE_PY" "$COVERAGE_SCRIPT" collect --manifest "$COVERAGE_MANIFEST" --output "$COVERAGE_RUN/collection_${SLURM_JOB_ID}.json"
  COVERAGE_AGREEMENT=$("$COVERAGE_PY" -c 'import json,sys;print(json.load(open(sys.argv[1]))["agreement"]["path"])' "$COVERAGE_MANIFEST")
  "$COVERAGE_PY" "$COVERAGE_SCRIPT" prepare_pool --proposals "$COVERAGE_RUN/collection_${SLURM_JOB_ID}.json" --agreement "$COVERAGE_AGREEMENT" --output "$COVERAGE_ROOT/pool_v1"
  "$COVERAGE_PY" "$COVERAGE_ROOT/pool_v1/implementation/preparation_coverage_score.py" execute_mace --manifest "$COVERAGE_ROOT/pool_v1/manifest.json"
else exit 2; fi
