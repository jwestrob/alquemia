# Six-source standalone solvent-response pilot

The actual frozen run is `workspaces/standalone_response_20260923/run_v3/`.
Job1210579 completed successfully once: all six pools, 228 standalone calls and
94 new MACE evaluations; no accuracy expansion justified. See [report](REPORT.md).
run_v2 is preparation-only and must not execute.
Historical four-source proposal inventory remains explicitly proposed-only;
PLAN_v2.md records the authorized six-source extension before any new output.

```bash
RESPONSE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
RESPONSE_RUN="$PWD/workspaces/standalone_response_20260923/run_v3"
"$RESPONSE_PY" "$RESPONSE_RUN/implementation/standalone_response_proposals.py" validate \
  --manifest "$RESPONSE_RUN/manifest.json"
```

The executed wrapper wrote the complete comparison automatically, preserving all
attempts. Do not submit it again to collect.
Its resource shape is32CPU/200000MiB/oneH200, with two endpoint searches and one
serialized warm native MACE worker. Both standalone media run as eight-thread jobs.

Report an actual completed collection (use a fresh output filename):

```bash
"$RESPONSE_PY" "$RESPONSE_RUN/implementation/standalone_response_proposals.py" report \
  --collection "$RESPONSE_RUN/collection.json" --output "$RESPONSE_RUN/review_comparison.json"
```

Raw low-level energy/gradient/state outputs are in baseline/, searches/ and cross/.
Each task has a manifest-bound execution receipt; every query and final proposal is
retained. STARTS.json records actual baseline-derived choices before optimization.
No native ORCA solvent scalar substitutes for a missing standalone cell.

Costs are already collected in COSTS.json. For a separate read-only audit, use fresh
accounting and output filenames:

```bash
sacct -j 1210579 --parsable2 --noheader --units=K \
  --format=JobIDRaw,State,ExitCode,AllocCPUS,ElapsedRaw,CPUTimeRAW,ReqMem,AllocTRES,MaxRSS \
  > "$RESPONSE_RUN/review_scheduler_accounting.txt"
"$RESPONSE_PY" diagnostics/standalone_response_20260923/summarize_execution.py \
  --run "$RESPONSE_RUN" --accounting "$RESPONSE_RUN/review_scheduler_accounting.txt" \
  --output "$RESPONSE_RUN/review_COSTS.json"
```

Actual-fixture tests (seven pass, zero skips after the completed pilot):

```bash
"$RESPONSE_PY" -m unittest discover -s tests -p test_standalone_response_proposals.py -v
```

No default, reference or deployment change; source selection is informed development.
