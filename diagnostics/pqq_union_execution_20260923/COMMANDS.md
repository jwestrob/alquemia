# Matched fresh execution: current preflight commands

The two molecular jobs are **not submitted**. Root will review the full225 result
or confirm the completed integration preflight before submission. Existing source
preparation/replay entrypoints and production remain unchanged.

From the repository root:

```bash
UNION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
UNION_BASE="$PWD/workspaces/pqq_union_execution_20260923"
"$UNION_PY" scripts/pqq_union_execution.py dry-run \
  --plan "$UNION_BASE/released-static_v4/plan.json"
"$UNION_PY" scripts/pqq_union_execution.py dry-run \
  --plan "$UNION_BASE/union-candidate_v4/plan.json"
"$UNION_PY" -m unittest discover -s tests -p test_pqq_union_execution.py -v
```

The source request is `workspaces/pqq_union_execution_20260923/REQUEST.json`:
all ten A0A3 sources, exact selectors and original canonical membership. Plans
v1/v2/v3 are retained unexecuted development snapshots; v4 is the ready snapshot.
`mapping_preflight_v4` uses real archived geometry solely to verify physical
mappings and generated solver inputs. It is outside the fresh execution paths,
contains no molecular outputs and cannot supply the fresh preparation/results.

When both approved jobs actually finish, their wrappers always collect/report,
including an execution failure. The matched report-only comparison is:

```bash
"$UNION_PY" "$UNION_BASE/union-candidate_v4/implementation/pqq_union_execution.py" compare \
  --static "$UNION_BASE/released-static_v4/RESULT.json" \
  --candidate "$UNION_BASE/union-candidate_v4/RESULT.json" \
  --output "$UNION_BASE/COMPARISON_v1.json"
```

No archived energy/force may fill a fresh result. The candidate can reuse its own
fresh q0 forces and optimizer final evaluations within the same run. Failed
preparation, origins, proposals or required cross-scored cells remain explicit.
The original static and candidate carry separate frozen references; the JSON
retains all ten rows and strict La4/Ca5/balanced/four La-triple summaries.

The batch wrapper is `diagnostics/pqq_union_execution_20260923/run.sbatch`.
Submit static first and candidate second, with the same permitted H200 host,
32 tasks/one CPU each, one GPU and200000MiB. A saved explicit submission receipt
must record the dependency, manifest/plan pins, executable wrapper, jobs and logs.
No time limit or scientific retry is added by this interface.
