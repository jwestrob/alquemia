# Direct donor-basin breadth operations

The grid and model are frozen in PLAN.md. These commands use explicit paths and
write-once products. Existing completed molecular jobs must not be resubmitted.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
BASIN_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
BASIN_WORK=workspaces/local_basin_breadth_20260922

# Existing real prepared grid replay; no molecular calls.
"$BASIN_PY" scripts/local_basin_breadth.py dry-run --design "$BASIN_WORK/grid_v1/design.json"

# Actual preparation commands already executed; retain these paths for provenance.
# prepare --inputs diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json \
#   --agreement diagnostics/local_basin_breadth_20260922/PLAN.md \
#   --output workspaces/local_basin_breadth_20260922/grid_v1
# specification --design "$BASIN_WORK/grid_v1/design.json" \
#   --output "$BASIN_WORK/SPECIFICATION_v1.json"

# Shared scoring adapter, not a second molecular execution system.
"$BASIN_PY" "$BASIN_WORK/scoring_v1/implementation/nikasha_finite_candidates.py" validate \
  --manifest "$BASIN_WORK/scoring_v1/manifest.json"

# Reuse the last existing complete collection; no repeated collection needed.
PYTHONPATH="$BASIN_WORK/scoring_v1/implementation" \
"$BASIN_PY" "$BASIN_WORK/analysis_implementation_v4/local_basin_breadth.py" analyze \
  --design "$BASIN_WORK/grid_v1/design.json" \
  --collection "$BASIN_WORK/scoring_v1/after_solvent_3_1210125.json" \
  --output "$BASIN_WORK/reanalyzed_v1.json"
PYTHONPATH="$BASIN_WORK/scoring_v1/implementation" \
"$BASIN_PY" "$BASIN_WORK/report_implementation_v1/local_basin_breadth.py" plot \
  --result "$BASIN_WORK/reanalyzed_v1.json" \
  --output "$BASIN_WORK/replotted_v1"
```

The scheduler submissions are captured as exact argument arrays in scoring_v1's
SUBMISSION_*.json receipts. Execution uses existing run_mace.sbatch and
run_solvent.sbatch wrappers, with four independent solvent manifests/locks.
All 65 nodes, including exact reused q0, remain required for each metal. The
analyzer excludes the non-grid base candidates retained by the scoring adapter.
It distinguishes raw conditional integrals from qualified numerical corrections.
No classifier calibration, whole-pocket entropy or chemical-state populations
are supplied by this operation.

Actual execution: MACE1210107 completed512/512 calls. Original standard-partition
allocations1210108–1210111 never started and were canceled after guarded checks;
replacement CPU-only GPU-partition jobs1210122–1210125 use the identical four
shard manifests. Full pre/post partition/node/cancellation/submission receipts
remain in scoring_v1. No scientific task was duplicated or retried.
