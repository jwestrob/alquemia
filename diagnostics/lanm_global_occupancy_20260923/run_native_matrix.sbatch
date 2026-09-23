#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J nikasha-lanm-solvent
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200000M
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
lanm_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
cd "$lanm_root"
lanm_python=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
lanm_execute_status=0
"$lanm_python" \
  workspaces/lanm_global_occupancy_20260923/scoring_v2/implementation/affordable_workflow.py execute \
  --manifest workspaces/lanm_global_occupancy_20260923/native_matrix_v2/manifest.json || lanm_execute_status=$?
lanm_collect_status=0
"$lanm_python" diagnostics/lanm_global_occupancy_20260923/continue_pilot.py collect-notify || lanm_collect_status=$?
"$lanm_python" - "$lanm_execute_status" "$lanm_collect_status" <<'PY'
import json,sys
from pathlib import Path
p=Path('workspaces/lanm_global_occupancy_20260923/NATIVE_PROCESS_STATUS.json')
with p.open('x') as f:
    json.dump({'executor_exit_status':int(sys.argv[1]),'collector_exit_status':int(sys.argv[2])},f,indent=2)
PY
if (( lanm_execute_status != 0 || lanm_collect_status != 0 )); then exit 1; fi
