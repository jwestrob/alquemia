#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-vacuum-esp
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
if [ "$#" -ne 3 ]; then
    echo 'Usage: run_esp.sbatch ABSOLUTE_MANIFEST ABSOLUTE_ESP_DIR ABSOLUTE_STATE_DIR' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$FIELD_PY" scripts/global_electrostatic.py esp --manifest "$1" --output "$2"
"$FIELD_PY" scripts/global_electrostatic.py prepare-states --manifest "$1" --esp-result "$2/esp_result.json" --output "$3"
