#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-xtb-reference
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
XTB_REFERENCE=$1
case "$XTB_REFERENCE" in /*) ;; *) exit 2;; esac
XTB_REF_RUN=$(dirname "$XTB_REFERENCE")
XTB_REF_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
ulimit -s unlimited
"$XTB_REF_PY" "$XTB_REF_RUN/implementation/standalone_xtb_reference.py" execute --manifest "$XTB_REFERENCE"
xtb_ref_status=$?
"$XTB_REF_PY" "$XTB_REF_RUN/implementation/standalone_xtb_reference.py" collect --manifest "$XTB_REFERENCE" --output "$XTB_REF_RUN/collection_v1.json"
exit "$xtb_ref_status"
