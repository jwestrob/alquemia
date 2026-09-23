# Primary225 proposal transfer: terminal result

All410 declared searches finished: **409 valid bounded candidates and one
explicit nonconverged search**. This gives204 complete endpoint pairs among205
eligible cases. The complete ledger retains all225 cases and450 endpoint rows,
including20 inherited unavailable cases; it does not substitute earlier scores
for failed candidates. Root owns the separately calibrated common-pool scoring
that determines classification performance.

## The one numerical failure

`q60ar6-pqq-la_model__conditioned_La__seed-1_sample-0__La` reached the unchanged
200-iteration SLSQP limit after2052 objective evaluations and201 gradient
evaluations. Search elapsed time was579.3921093419194 seconds. Native model
receipts continued to advance, and the worker proceeded to subsequent endpoints
normally. This was optimizer nonconvergence, not a stalled GPU or missing source
chemistry. The candidate remains unavailable; no extension, alternate scale,
restart, relaxed threshold or earlier-score fallback was applied.

Among409 valid candidates,172 reach an angular or displacement boundary. They
are constrained geometry proposals, not demonstrated unconstrained minima,
ensemble populations or free energies. The same original-q0 mode selector,
physical atoms, charges, water inventory, one-start policy and native-MACE
energy/force model were applied to every eligible case.

All accepted native energy changes are nonpositive. Maximum accepted heavy-atom
displacement is0.8000000000228809 Å, within the unchanged1e−7 Å tolerance.
Intermediate trials generated1199 completed requests outside the final0.8 Å
domain, with maximum2.0595045512210826 Å; repeated requests may reuse a point.
As declared, the displacement limit applies to final candidates while unchanged
chemical/overlap guards apply throughout.

The searches produced11,541 actual new native-MACE result receipts. Median
per-search elapsed time is4.751154217869043 seconds; the95th percentile is
13.610064441151712 seconds. The579.392-second nonconverged search is a substantial
outlier. These search intervals exclude shared allocation/validation overhead
and future composite scoring; they are not complete per-site production costs.

## Actual allocation cost

| Shard/job | Searches | Allocation seconds | Allocated CPU-s | Allocated GPU-s |
|---|---:|---:|---:|---:|
| 0 /1210021 | 103 | 735 | 23,520 | 735 |
| 1 /1210022 | 103 | 807 | 25,824 | 807 |
| 2 /1210023 | 102 | 606 | 19,392 | 606 |
| 3 /1210024 | 102 | 1,311 | 41,952 | 1,311 |
| Total | **410** | **3,459 summed** | **110,688** | **3,459** |

Each allocation reserved32 CPUs,1 H200 and200000 MiB. All four Slurm jobs
completed successfully; a job can complete while preserving an individual
scientific task's explicit failure. Actual elapsed time from first allocation
start to final end was1,940 seconds (32 min20 s), with peak concurrency3 GPUs.
Allocation sums above include the failed search and worker overhead. They do
not misrepresent overlap as lower total compute cost. No DFT/GFN2 calls were
made by these jobs. Model q0 values/forces were reused from exact source receipts.

## Artifacts and interpretation

All products are under
`workspaces/adaptive_completion_20260922/primary225_v2_sharded/`:

- `manifest.json`: immutable410-task source/optimizer pins, SHA256
  `6d94824785d630d294f751b36aa86edbb90ab4cdb7962f8a68658344e5e7ba1c`.
- `SUBMISSION.json`, four `execution_JOB_shard_N.json` receipts and
  `ACCOUNTING.txt`: exact execution selection and actual scheduler accounting.
- `final_collection.json`: full225/450 scientific ledger after physical mapping
  validation; SHA256
  `5c674a8674951f213229bfef82ec2a01928b5eab18b418f3d8b321e9374d14c0`.
  The incomplete pair and inherited missing inputs remain explicit.
- `TERMINAL_SUMMARY.json` and `summarize_terminal.py`: all410 terminal rows,
  workload counts, timing quantiles and exact source/receipt/summary-code pins.
- Per-task `result.json`, Cartesian forces, model receipts and optimizer traces
  preserve all evaluated trials and terminal acceptance decisions.

The existing canonical-only reference was frozen before transfer execution.
These noncanonical folds were previously inspected development data, so this
is structural transfer evidence, not prospectively blind validation. Classification
benefit and total scoring affordability require root's subsequent cross-metal
and solvent calculation. Baseline/production remain unchanged.
