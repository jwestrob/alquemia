#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J nikasha-lanm-native-retry
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=200000M
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
lanm_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
cd "$lanm_root"
lanm_py=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
lanm_exit=0
"$lanm_py" workspaces/lanm_global_occupancy_20260923/scoring_v2/implementation/affordable_workflow.py execute \
 --manifest workspaces/lanm_global_occupancy_20260923/native_feasibility_retry_v1/manifest.json || lanm_exit=$?
"$lanm_py" diagnostics/lanm_global_occupancy_20260923/collect_native_recovery.py
exit "$lanm_exit"
