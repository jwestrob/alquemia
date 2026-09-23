# Opt-in fresh union/adaptive execution

The matched fresh trial is complete: [report](REPORT.md), [pins](ARTIFACTS.json).
Both results already exist; reading them requires no molecular work. Production
remains unchanged. This narrow execution adapter currently accepts the exact
all-ten A0A3F2YLY8 reference-source request; it is not a generic unknown-protein or
three-fold interface. The older archive-replay entrypoint remains unchanged.

## Read the completed result

From the repository root:

```bash
cat diagnostics/pqq_union_execution_20260923/REPORT.md
cat workspaces/pqq_union_execution_20260923/COMPARISON_v1.json
```

The exact original request is
`workspaces/pqq_union_execution_20260923/REQUEST.json`: all five Ca-conditioned
and five La-conditioned source files, their selectors/configuration and actual
canonical membership. Final collections are `{released-static_v4,union-candidate_v4}/RESULT.json`.
The static collection-only recovery is explicitly pinned in
`released-static_v4/COLLECTION_RECOVERY_v2.json`. V1–v3 plans remain unexecuted;
v4 scientific snapshots remain unchanged after execution.

## Prepare and execute an explicit fresh repeat

These commands create a **new** sandbox and perform fresh work if submitted.
No archived energy/force fills a cell. The candidate's maximum is20 fresh origin
MACE evaluations,20 bounded searches,20 cross-MACE evaluations and120 GFN2 cells.
Use the completed results above when another calculation is unnecessary.

```bash
UNION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
UNION_REPEAT="$PWD/workspaces/pqq_union_execution_repeat_v1"
"$UNION_PY" scripts/pqq_union_execution.py prepare \
  --request workspaces/pqq_union_execution_20260923/REQUEST.json \
  --arm union-candidate \
  --reference workspaces/slsqp_precision_expansion_20260923/REFERENCE_v1.json \
  --rank-qualification workspaces/native_gfn2_rank_panel_20260923/run_v2/COMPARISON.json \
  --maxiter-qualification workspaces/adaptive_maxiter_20260922/prepared_v1/collection_1209876.json \
  --agreement diagnostics/pqq_union_execution_20260923/EXECUTION_AGREEMENT.md \
  --output "$UNION_REPEAT"
"$UNION_PY" "$UNION_REPEAT/implementation/pqq_union_execution.py" dry-run \
  --plan "$UNION_REPEAT/plan.json"
sbatch --parsable --partition=gpu --nodelist=node-224-2t-8gpu-1 \
  --nodes=1 --ntasks=32 --cpus-per-task=1 --gres=gpu:1 --mem=200000M \
  --job-name=nikasha-optin-union \
  --output "$UNION_REPEAT/slurm_%j.out" --error "$UNION_REPEAT/slurm_%j.err" \
  diagnostics/pqq_union_execution_20260923/run.sbatch "$UNION_REPEAT/plan.json"
```

The wrapper runs the pinned implementation, then always collects and reports,
even after an execution failure. Save the actual scheduler receipt when using
this manual operation. It writes `RESULT.json` and `RESULT.md`, or job-specific
names if a previous result exists. Manual report-only operations, using a new
output name to preserve prior records, are:

```bash
"$UNION_PY" scripts/pqq_union_execution.py collect \
  --plan "$UNION_REPEAT/plan.json" --output "$UNION_REPEAT/RESULT_manual_v1.json"
"$UNION_PY" scripts/pqq_union_execution.py report \
  --result "$UNION_REPEAT/RESULT_manual_v1.json" --output "$UNION_REPEAT/REPORT_manual_v1.md"
```

The live collector includes the narrow top-level scanner-filename fix discovered
in the completed trial; no executed input/output or snapshot was modified. The
candidate retains all ten rows, original union R, accommodated R, all matrix
components, selections/works, boundaries, states and missing reasons. The separate
released static comparison uses `--arm released-static` with its own new output
directory, frozen released bands and original solver path. Exact v4 submissions
for both sequential arms are saved under their `SUBMISSION.json` records.

## Compare, cost and verify the completed trial

The following operations are report-only; output paths are new versions:

```bash
"$UNION_PY" scripts/pqq_union_execution.py compare \
  --static workspaces/pqq_union_execution_20260923/released-static_v4/RESULT.json \
  --candidate workspaces/pqq_union_execution_20260923/union-candidate_v4/RESULT.json \
  --output workspaces/pqq_union_execution_20260923/COMPARISON_manual_v1.json
"$UNION_PY" -m unittest discover -s tests -p test_pqq_union_execution.py -v
```

Actual scalar/native costs and source/score reproduction are retained in
`COSTS_v2.json` and `REPRODUCTION_v1.json`. Their explicit-path report helpers are
`summarize_costs.py` and `compare_archives.py` in this directory. Whole allocation
accounting includes all stages; stage timings must not be double-counted.
[Timing interpretation](TIMING_NOTE.md), [restart limits](RESTART_NOTE.md).

Failed preparation or required origin/proposal/matrix cells remain unavailable.
No failed molecule is retried automatically. A started released static runner is
collection-only on restart; it has no stage-resume API. These limits do not change
the successful scientific results or authorize an alternate source/threshold.
