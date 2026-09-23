# Restart and collective scaffold operations

Run from the repository root. The saved submissions have already been executed;
do not resubmit them. Commands below read artifacts or create fresh reports.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
NIKASHA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
NIKASHA_POOL=$PWD/workspaces/collective_scaffold_20260922/pool_v1
```

Audit the exact scored candidate preparation without molecular calculations:

```bash
"$NIKASHA_PY" "$NIKASHA_POOL/implementation/nikasha_finite_candidates.py" validate \
  --manifest "$NIKASHA_POOL/manifest.json"
```

Repeat the final comparison into a new directory, leaving published results intact:

```bash
NIKASHA_REPORT=$(mktemp -d workspaces/scaffold-comparison.XXXXXX)
"$NIKASHA_PY" "$NIKASHA_POOL/implementation/nikasha_parallel_compare.py" \
  --inputs diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json \
  --collections "$NIKASHA_POOL/final_collection.json" \
  --output "$NIKASHA_REPORT/COMPARISON.json"
```

Actual candidate preparation is pinned in `SCORING_SPEC_v1.json` alongside the
pool directory. `nikasha_scaffold_candidates.py` creates that specification from
one complete mechanics collection using explicit `--collection`, `--inputs`,
`--agreement` and `--output` arguments. The existing finite-candidate `prepare`
operation takes `--specification`, `--output` and `--shards 1`. It reconstructs and
checks all physical candidates, writes the existing worker manifests and snapshots
the scientific implementation. No new source paths are hidden in the scorer.

The prepared pool's `PREFLIGHT.json` is the dry-run record. Exact contained execute
commands, allocations and successful job IDs are in `SUBMISSION_MACE.json` and
`SUBMISSION_solvent.json`. These use the existing common-pool wrappers. The CPU
solvent submission explicitly selected the gpu partition's CPU-sharing route,
with64 allocated CPUs/eight8-rank workers and no GPU request. MACE uses oneH200,
32CPUs and200000MiB. Neither wrapper imposes a new project-total time/compute cap.
Each collects its existing partial results; final collection requires both jobs.

To recollect the existing endpoint outputs into the fresh report directory:

```bash
"$NIKASHA_PY" "$NIKASHA_POOL/implementation/nikasha_pool.py" collect \
  --manifest "$NIKASHA_POOL/manifest.json" \
  --output "$NIKASHA_REPORT/recollected.json"
```

A failed required cell remains unavailable. Native continuation
commands and its actual first/recovery receipts are documented separately in
`diagnostics/native_xtb_restart_20260922/COMMANDS.md`. Scaffold mechanics version3
commands and all retained failed/superseded attempts are in
`diagnostics/collective_scaffold_20260922/COMMANDS.md`.

The final comparison reports every required matched target, q0 response control,
expanded pool, both transferred historical band sets, two prescribed structural
pair spreads and actual molecular call counts. It neither refits a threshold nor
promotes the candidate. Parent-force diagnostics and quantum restart diagnostics
remain separate from a full La/Ca scoring claim.
