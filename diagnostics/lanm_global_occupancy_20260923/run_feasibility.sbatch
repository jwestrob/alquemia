#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-lanm-global
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
lanm_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
cd "$lanm_root"
lanm_impl="$lanm_root/workspaces/lanm_global_occupancy_20260923/scoring_v2/implementation"
lanm_workspace="$lanm_root/workspaces/lanm_global_occupancy_20260923"
# At most 24 MPI ranks plus 8 MACE host threads within this 32-CPU allocation.
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  "$lanm_impl/affordable_workflow.py" execute \
  --manifest "$lanm_workspace/native_feasibility_v2/manifest.json" \
  > "$lanm_workspace/native_feasibility_v2/executor.log" 2>&1 &
lanm_native_pid=$!
lanm_gpu_status=0
"$lanm_root/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" \
  "$lanm_impl/lanm_global_occupancy.py" mace \
  --manifest "$lanm_workspace/scoring_v2/manifest.json" \
  --state Hans_8DQ2__EF12 --output "$lanm_workspace/mace_feasibility_v2" || lanm_gpu_status=$?
lanm_native_status=0
wait "$lanm_native_pid" || lanm_native_status=$?
if (( lanm_gpu_status != 0 || lanm_native_status != 0 )); then
  exit 1
fi
