# Native GFN2 rank experiment

All 24 cells completed in immutable `run_v1`: jobs1210743/1210744/1210745
for1/4/8ranks, sequential afterany. The saved SUBMISSION_rank_N.json records
each explicit CPU-only allocation and command. Do not resubmit to collect.

```bash
RANK_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
RANK_RUN="$PWD/workspaces/native_gfn2_ranks_20260923/run_v1"
"$RANK_PY" "$RANK_RUN/implementation/native_gfn2_rank_scaling.py" validate \
  --manifest "$RANK_RUN/rank_1/manifest.json"
"$RANK_PY" -m unittest discover -s tests -p test_native_gfn2_rank_scaling.py -v
```

Each wrapper collects its real task statuses even when execution fails. The
eight-rank wrapper then writes COMPARISON.json using all three collections and
explicit missing values. It does not repeat any molecular call.

Read-only comparison rerun, with a fresh output filename:

```bash
"$RANK_PY" "$RANK_RUN/implementation/native_gfn2_rank_scaling.py" compare \
  --design "$RANK_RUN/design.json" --output "$RANK_RUN/review_COMPARISON.json"
```

The terminal accounting and COSTS.json are saved. A fresh read-only accounting
and cost audit can use different output filenames:

```bash
sacct -j 1210743,1210744,1210745 --parsable2 --noheader --units=K \
  --format=JobIDRaw,State,ExitCode,AllocCPUS,ElapsedRaw,CPUTimeRAW,TotalCPU,ReqMem,AllocTRES,MaxRSS \
  > "$RANK_RUN/review_scheduler_accounting.txt"
"$RANK_PY" diagnostics/native_gfn2_ranks_20260923/summarize_costs.py \
  --design "$RANK_RUN/design.json" --accounting "$RANK_RUN/review_scheduler_accounting.txt" \
  --output "$RANK_RUN/review_COSTS.json"
```

The dataset and electronic inputs are fixed by [PLAN.md](PLAN.md), the design,
source manifest and exact archive receipts. Rank count is an execution
diagnostic; no result or shared runtime is promoted by these commands.
