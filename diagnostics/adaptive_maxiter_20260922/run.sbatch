#!/bin/bash
#SBATCH -p gpu
#SBATCH --nodelist=node-224-2t-8gpu-1
#SBATCH -J nikasha-maxiter-check
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
set -uo pipefail
DIAG_MANIFEST=$1
case "$DIAG_MANIFEST" in /*) ;; *) exit 2;; esac
DIAG_RUN=$(dirname "$DIAG_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
DIAG_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$DIAG_PY" "$DIAG_RUN/implementation/adaptive_maxiter_diagnostic.py" execute --manifest "$DIAG_MANIFEST"
DIAG_EXIT=$?
"$DIAG_PY" "$DIAG_RUN/implementation/adaptive_maxiter_diagnostic.py" collect --manifest "$DIAG_MANIFEST" --output "$DIAG_RUN/collection_${SLURM_JOB_ID}.json"
COLLECT_EXIT=$?
if [ "$DIAG_EXIT" -ne 0 ]; then exit "$DIAG_EXIT"; fi
exit "$COLLECT_EXIT"
