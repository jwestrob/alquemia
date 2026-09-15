#!/bin/bash
#SBATCH --partition=high-memory
#SBATCH --job-name=pqq_1kb0_fc_v3
#SBATCH --output=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_q46444_1kb0_external_validation_20260915/slurm_external_%j.out
#SBATCH --error=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_q46444_1kb0_external_validation_20260915/slurm_external_%j.err

set -euo pipefail

HERE=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_q46444_1kb0_external_validation_20260915
PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
ORCA_DIR=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
PREPARATION=${1:-"$HERE/prepared/external_preparation.json"}

export ORCA_PATH="$ORCA_DIR"
export PATH="$ORCA_DIR:$PATH"
export LD_LIBRARY_PATH="$ORCA_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$HERE"
"$PYTHON" "$HERE/run_external.py" \
  --preparation "$PREPARATION" \
  --pins "$HERE/implementation_pins.json"

