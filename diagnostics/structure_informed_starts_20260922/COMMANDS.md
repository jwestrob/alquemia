# Reproduce the completed alternative-start pilot

All four scientific jobs (1210098, 1210103, 1210104, 1210115) are complete. Do not
submit these experiments again. Current results and the follow-on proposal are separate.
Commands below are read-only validation or collection from existing outputs.
Every new destination must not already exist.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
STARTS_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
STARTS_ROOT=$PWD/workspaces/structure_informed_starts_20260922
cat diagnostics/structure_informed_starts_20260922/RESULT.json
"$STARTS_PY" "$STARTS_ROOT/searches_v1/implementation/structure_informed_starts.py" validate --manifest "$STARTS_ROOT/searches_v1/manifest.json"
"$STARTS_PY" "$STARTS_ROOT/pool_v1/implementation/nikasha_finite_candidates.py" validate --manifest "$STARTS_ROOT/pool_v1/manifest.json"
"$STARTS_PY" "$STARTS_ROOT/numerical_v2/implementation/affordable_workflow.py" dry-run --manifest "$STARTS_ROOT/numerical_v2/manifest.json"
```

Recollect the actual molecular outputs into a new artifact:

```bash
"$STARTS_PY" "$STARTS_ROOT/searches_v1/implementation/structure_informed_starts.py" collect --manifest "$STARTS_ROOT/searches_v1/manifest.json" --output "$STARTS_ROOT/searches_v1/collection_replay.json"
"$STARTS_PY" "$STARTS_ROOT/pool_v1/implementation/nikasha_pool.py" collect --manifest "$STARTS_ROOT/pool_v1/manifest.json" --output "$STARTS_ROOT/pool_v1/collection_replay.json"
"$STARTS_PY" "$STARTS_ROOT/numerical_v2/implementation/structure_informed_starts.py" collect-numeric --manifest "$STARTS_ROOT/numerical_v2/manifest.json" --output "$STARTS_ROOT/numerical_v2/collection_replay.json" > "$STARTS_ROOT/numerical_v2/collection_replay.log"
PYTHONPATH="$PWD/scripts" "$STARTS_PY" "$STARTS_ROOT/analysis_v2/structure_informed_starts.py" compare --collection "$STARTS_ROOT/pool_v1/final_collection.json" --searches "$STARTS_ROOT/searches_v1/collection_1210098.json" --output "$STARTS_ROOT/analysis_v2/comparison_replay.json" > "$STARTS_ROOT/analysis_v2/comparison_replay.log"
```

The comparison uses immutable existing static/adaptive bands and preserves the
Q9 edge crossing exactly. It creates no new classifier calibration. Missing starts
remain unavailable, even though their earlier scores still exist separately.

Geometry-only recreation and actual-fixture tests:

```bash
"$STARTS_PY" scripts/structure_informed_starts.py prepare --inputs diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json --agreement diagnostics/structure_informed_starts_20260922/PLAN.md --output "$STARTS_ROOT/preparation_replay"
"$STARTS_PY" -m unittest discover -s tests -p test_structure_informed_starts.py -v
```

The successful scientific submissions are recorded exactly in:

- searches_v1/SUBMISSION.json (wrapper run_proposals.sbatch).
- pool_v1/SUBMISSION_mace.json and SUBMISSION_solvent.json (existing pool wrappers).
- numerical_v2/SUBMISSION.json (wrapper run_numerical.sbatch).

The unsubmitted numerical_v1/Tight proposal is historical and must not be launched
as a qualified precision test. See NUMERICAL_REPEAT_PLAN.md for the documented
native-mixer limitation. RESTART_PROPOSED.md and RESTART_SOURCES.json describe a
possible separate native .xtbw continuity test; they are not an execution manifest.
