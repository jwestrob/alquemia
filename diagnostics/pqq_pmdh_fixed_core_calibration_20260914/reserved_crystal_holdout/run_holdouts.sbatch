#!/bin/bash
#SBATCH --partition=high-memory
#SBATCH --job-name=pqq_fc_holdout_v3
#SBATCH --output=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/slurm_holdout_%j.out
#SBATCH --error=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/slurm_holdout_%j.err

set -euo pipefail

HERE=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout
PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
ORCA_DIR=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
PREPARATION=${1:?usage: sbatch run_holdouts.sbatch /absolute/path/holdout_preparation.json [/fresh/result/dir]}
RECEIPT="$(dirname -- "$PREPARATION")/holdout_execution_${SLURM_JOB_ID}.json"
RESULT_DIR=${2:-"$HERE/result"}

export ORCA_PATH="$ORCA_DIR"
export PATH="$ORCA_DIR:$PATH"
export LD_LIBRARY_PATH="$ORCA_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$HERE"
if [[ -e "$RESULT_DIR" ]]; then
  echo "ERROR: refusing to overwrite existing result path: $RESULT_DIR" >&2
  exit 2
fi

set +e
"$PYTHON" "$HERE/run_holdouts.py" \
  --preparation "$PREPARATION" \
  --holdout-pins "$HERE/holdout_implementation_pins.json"
RUNNER_STATUS=$?
set -e
echo "holdout runner exit status: $RUNNER_STATUS"

if [[ ! -f "$RECEIPT" ]]; then
  echo "ERROR: runner exited $RUNNER_STATUS without writing receipt: $RECEIPT" >&2
  exit 2
fi

"$PYTHON" "$HERE/score_holdouts.py" \
  --preparation "$PREPARATION" \
  --execution-receipt "$RECEIPT" \
  --output-dir "$RESULT_DIR" \
  --holdout-pins "$HERE/holdout_implementation_pins.json"
