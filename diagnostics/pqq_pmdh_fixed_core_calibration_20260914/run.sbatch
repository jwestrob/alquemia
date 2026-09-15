#!/bin/bash
#SBATCH -p high-memory
#SBATCH -J pqq_fc_cal_v3
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_pmdh_fixed_core_calibration_20260914/slurm_%j.out
#SBATCH -e /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_pmdh_fixed_core_calibration_20260914/slurm_%j.err
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14

set -euo pipefail

HERE=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_pmdh_fixed_core_calibration_20260914
PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
ORCA_DIR=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg

export ORCA_PATH="$ORCA_DIR"
export PATH="$ORCA_DIR:$PATH"
export LD_LIBRARY_PATH="$ORCA_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd "$HERE"
"$PYTHON" "$HERE/run_panel.py" \
  --preparation "$HERE/prepared/preparation.json" \
  --implementation-pins "$HERE/implementation_pins.json"

"$PYTHON" "$HERE/score.py" \
  --preparation "$HERE/prepared/preparation.json" \
  --implementation-pins "$HERE/implementation_pins.json" \
  --output-dir "$HERE"
