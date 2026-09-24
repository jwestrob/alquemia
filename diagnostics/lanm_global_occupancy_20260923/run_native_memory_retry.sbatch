#!/bin/bash
#SBATCH -p memory
#SBATCH -J nikasha-lanm-native-memory
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jacobwestroberts@gmail.com
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
lanm_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
cd "$lanm_root"
lanm_py=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
lanm_output=workspaces/lanm_global_occupancy_20260923/native_feasibility_retry_v2
# Sixty-four is the scheduler's minimum task request; exclusive grants the whole
# node. The runtime manifest uses ALL allocated CPUs, including on larger nodes.
read -r lanm_ranks lanm_memory_mib lanm_maxcore_mb < <("$lanm_py" diagnostics/lanm_global_occupancy_20260923/native_allocation.py)
"$lanm_py" diagnostics/lanm_global_occupancy_20260923/prepare_native_memory_retry.py \
 --source-manifest workspaces/lanm_global_occupancy_20260923/native_feasibility_retry_v1/manifest.json \
 --output "$lanm_output" --ranks "$lanm_ranks" --memory-mib "$lanm_memory_mib" --maxcore-mb "$lanm_maxcore_mb"
lanm_manifest="$lanm_output/manifest.json"
lanm_exit=0
"$lanm_py" workspaces/lanm_global_occupancy_20260923/scoring_v2/implementation/affordable_workflow.py execute \
 --manifest "$lanm_manifest" || lanm_exit=$?
"$lanm_py" diagnostics/lanm_global_occupancy_20260923/collect_native_recovery.py --manifest "$lanm_manifest"
exit "$lanm_exit"
