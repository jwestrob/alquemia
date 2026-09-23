#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-standalone-xtb
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
XTB_MANIFEST=$1
case "$XTB_MANIFEST" in /*) ;; *) exit 2;; esac
XTB_RUN=$(dirname "$XTB_MANIFEST")
XTB_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
ulimit -s unlimited
"$XTB_PY" "$XTB_RUN/implementation/standalone_xtb.py" execute --manifest "$XTB_MANIFEST"
xtb_execution_status=$?
"$XTB_PY" "$XTB_RUN/implementation/standalone_xtb.py" collect --manifest "$XTB_MANIFEST" --output "$XTB_RUN/collection_v1.json"
exit "$xtb_execution_status"
