#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-xtb-static225
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
XTB_TRANSFER=$1
case "$XTB_TRANSFER" in /*) ;; *) exit 2;; esac
XTB_TRANS_RUN=$(dirname "$XTB_TRANSFER")
XTB_TRANS_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
ulimit -s unlimited
"$XTB_TRANS_PY" "$XTB_TRANS_RUN/implementation/standalone_xtb_transfer.py" execute --manifest "$XTB_TRANSFER"
xtb_trans_status=$?
"$XTB_TRANS_PY" "$XTB_TRANS_RUN/implementation/standalone_xtb_transfer.py" collect --manifest "$XTB_TRANSFER" --output "$XTB_TRANS_RUN/collection_v1.json"
exit "$xtb_trans_status"
