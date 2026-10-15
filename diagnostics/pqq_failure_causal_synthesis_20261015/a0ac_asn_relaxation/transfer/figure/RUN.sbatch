#!/bin/bash
#SBATCH --job-name=Asn-transfer-figure
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
if [ "$#" -lt 4 ]; then
  echo 'Usage: RUN.sbatch ABS_PLOT_PY SCRIPT_SHA ABS_REPORT_DIR ABS_NEW_OUTPUT_DIR [operational|mathematical]' >&2
  exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11 "$1" \
  --script-sha "$2" --report-dir "$3" --output "$4" --policy "${5:-operational}"
