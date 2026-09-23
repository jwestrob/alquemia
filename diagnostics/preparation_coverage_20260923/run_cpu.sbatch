#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-H-repair-GFN
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --nodelist=node-128-512g-8gpu-1
set -uo pipefail
COVERAGE_MANIFEST=$1
COVERAGE_PHASE=$2
case "$COVERAGE_MANIFEST" in /*) ;; *) exit 2;; esac
COVERAGE_RUN=$(dirname "$COVERAGE_MANIFEST")
COVERAGE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$COVERAGE_PY" "$COVERAGE_RUN/implementation/affordable_workflow.py" execute --manifest "$COVERAGE_RUN/solvent/shard_0/manifest.json"
COVERAGE_EXEC=$?
if [ "$COVERAGE_PHASE" = origins ]; then COVERAGE_OP=collect_origins; elif [ "$COVERAGE_PHASE" = pool ]; then COVERAGE_OP=collect_pool; else exit 2; fi
"$COVERAGE_PY" "$COVERAGE_RUN/implementation/preparation_coverage_score.py" "$COVERAGE_OP" --manifest "$COVERAGE_MANIFEST" --output "$COVERAGE_RUN/collection_final.json"
COVERAGE_COLLECT=$?
if [ "$COVERAGE_EXEC" -ne 0 ]; then exit "$COVERAGE_EXEC"; fi
exit "$COVERAGE_COLLECT"
