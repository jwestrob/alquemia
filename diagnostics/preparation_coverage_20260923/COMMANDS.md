# Supplemental repaired source: exact operations

No production default or historical225 row changes. The original16 misplaced
ions remain unsupported. Use immutable manifests/snapshots below; do not submit
jobs already recorded in SUBMISSION files.

```bash
coverage_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
coverage_python=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
coverage_run="$coverage_root/workspaces/preparation_coverage_20260923"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$coverage_python" "$coverage_run/origins_v1/implementation/preparation_coverage_score.py" validate_origins \
  --manifest "$coverage_run/origins_v1/manifest.json"
```

GPU1210669 owns four origins. Dependent CPU1210670 owns the eight native GFN2
cells and writes origins_v1/collection_final.json. Both use the unchanged finite
source manifest SHA41fd83962bff17c7407f8477fb36c4b782ab4b0c3cafe8fecfe0d8ca2fe0e43a.

After those actual outputs complete, the source snapshot constructs the two
bounded searches from the real union forces and the unchanged original policy:

```bash
"$coverage_python" "$coverage_run/origins_v1/implementation/preparation_coverage_score.py" prepare_proposals \
  --origins "$coverage_run/origins_v1/collection_final.json" \
  --output "$coverage_run/proposals_v1"
"$coverage_python" "$coverage_run/proposals_v1/implementation/preparation_coverage_score.py" validate \
  --manifest "$coverage_run/proposals_v1/manifest.json"
```

The GPU wrapper's `proposals` phase executes exactly two starts, collects actual
candidates, prepares the common pool and executes only missing cross-MACE cells.
The CPU wrapper's `pool` phase executes the actual finite solvent manifest and
collects it. Submit only once; record both submission receipts. A failed q0,
search or molecular cell remains explicit and never becomes a baseline fallback.

Report-only comparison after the actual pool exists:

```bash
"$coverage_python" "$coverage_root/scripts/preparation_coverage_score.py" compare \
  --origins "$coverage_run/origins_v1/collection_final.json" \
  --pool "$coverage_run/pool_v1/collection_final.json" \
  --local-reference "$coverage_root/diagnostics/pqq_fast_release_20260920/RELEASE_PANEL.json" \
  --union-reference "$coverage_root/diagnostics/consistent_context_20260922/FROZEN_REFERENCE_v1.json" \
  --output "$coverage_run/COMPARISON_v1.json"
"$coverage_python" -m unittest discover -s "$coverage_root/tests" -p test_preparation_coverage.py -v
```

Comparison records its actual implementation hash; the molecular source snapshot
predates the addition of this reporting helper. No scientific runtime changes.
The released artifact's frozen_bands are composite reference bands. The older
full_v1 comparison's top-level native_PQQ_bands describe vacuum-only results;
its calibration.context.bands hold the composite calibration. Do not confuse
those quantities or silently recalculate the released threshold.

## Completed delivery

All four allocations are terminal:1210669,1210670,1210819,1210820. Do not rerun
the molecular commands above. All16GFN and42 actualMACE calls succeeded. Read
REPORT.md/RESULT.json and the pinned COMPARISON_v1.json for the final outcome;
comparison outputs are write-once. The read-only fixture command above now
checks all completed actual artifacts (7tests,0skips).
