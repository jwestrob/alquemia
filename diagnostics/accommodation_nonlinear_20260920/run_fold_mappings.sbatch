#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J primary-fold-donor-maps
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
MAP_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
MAP_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$MAP_PY" "$MAP_ROOT/workspaces/accommodation_nonlinear_20260920/fold_maps_implementation_v1/accommodation_fold_mappings.py" \
  --preparation "$MAP_ROOT/workspaces/accommodation_goal_20260920/folds_v1/preparation_reconciled_v1.json" \
  --existing "$MAP_ROOT/workspaces/accommodation_reference_geometry_20260920/mappings_v1/result.json" \
  --engine-manifest "$MAP_ROOT/workspaces/accommodation_goal_20260920/folds_v1/solvent_shards_v2/shard_0/manifest.json" \
  --agreement "$MAP_ROOT/diagnostics/accommodation_nonlinear_20260920/FOLD_TRANSFER_PLAN.md" \
  --output "$MAP_ROOT/workspaces/accommodation_nonlinear_20260920/fold_maps_v1" --workers 64
