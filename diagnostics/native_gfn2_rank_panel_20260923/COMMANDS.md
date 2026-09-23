# Fixed 28-source rank qualification

Run from the repository root. These operations do not change any released default.
The actual run directory is `workspaces/native_gfn2_rank_panel_20260923/run_v2`.
`run_v1` is an unexecuted preparation attempt retained after correcting support
for the older native worker's successful `computed` status.

```bash
PANEL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
PANEL_ROOT="$PWD/workspaces/native_gfn2_rank_panel_20260923/run_v2"
PANEL_SCRIPT="$PANEL_ROOT/implementation/native_gfn2_rank_panel.py"

"$PANEL_PY" "$PANEL_SCRIPT" validate --manifest "$PANEL_ROOT/manifest.json"

# Submission is recorded in run_v2/SUBMISSION.json; do not submit a duplicate.
# The batch wrapper always collects completed and failed cells after execution.

# Read the immutable scientific verdict once the batch has finished:
"$PANEL_PY" -c 'import json,sys; r=json.load(open(sys.argv[1])); print({k:r[k] for k in ("complete","case_denominator","complete_cells","cell_denominator","all_numerical_gates_pass")})' "$PANEL_ROOT/COMPARISON.json"

"$PANEL_PY" -m unittest discover -s tests -p test_native_gfn2_rank_panel.py -v
```

The preparation command (already run) was:

```bash
"$PANEL_PY" scripts/native_gfn2_rank_panel.py prepare \
  --reference workspaces/union_adaptive_20260923/CANONICAL_REFERENCE_v1.json \
  --pilot workspaces/native_gfn2_ranks_20260923/run_v1/rank_1/collection.json \
  --agreement diagnostics/native_gfn2_rank_panel_20260923/PLAN.md \
  --output workspaces/native_gfn2_rank_panel_20260923/run_v2
```

For an interrupted *collection only*, first check for a live batch/collector and
existing output. `collect --manifest … --output …` and
`compare --collection … --output …` require fresh output paths and never launch
molecular calculations. Execution itself has no automatic failed-cell rescue.
