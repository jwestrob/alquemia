#!/bin/bash
#SBATCH -p standard
#SBATCH -J nikasha-union-solvent
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
set -euo pipefail
UNION_MANIFEST=$1
case "$UNION_MANIFEST" in /*) ;; *) exit 2;; esac
UNION_ROOT=$(dirname "$(dirname "$UNION_MANIFEST")")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNION_CPU_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$UNION_CPU_PY" "$UNION_ROOT/implementation/consistent_context.py" validate-solvent --manifest "$UNION_MANIFEST"
"$UNION_CPU_PY" "$UNION_ROOT/implementation/affordable_workflow.py" execute --manifest "$UNION_MANIFEST"
