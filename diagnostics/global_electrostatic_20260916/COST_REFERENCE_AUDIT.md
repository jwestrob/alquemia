# Cost reference audit: matching 1H4I cores

**2026-09-16. Exact no-MBIS baseline receipts and CPCM+MBIS pilot receipts
exist for both qm33/qm36 La/Ca pairs.** All four XYZ hashes match the new vacuum
inputs. ORCA is the same 6.1.1 executable. Hardware, MPI sizes and concurrency
are not fully matched, so these records do **not** establish a controlled
MBIS-overhead ratio or a production-cost ratio.

The main concrete finding is in ORCA's own timings: population analysis takes
**542.740–618.369 s** in the old CPCM+MBIS endpoints and **574.688–747.686 s**
in the new vacuum+MBIS endpoints. The archived no-MBIS runs report
**0.371–1.740 s** for that category. Population analysis includes more than
MBIS alone, so these are measured category times, not an isolated MBIS timer.
They locate a substantial cost in charge analysis. ORCA includes this work in
its larger SCF module; the final “SCF iterations” percentage must not be read
as electronic convergence alone.

## Endpoint receipts

Wall time below is finish minus start from the execution receipt. Rank-wall
seconds are that interval times configured MPI ranks (one thread/rank), **not
measured CPU usage or a separately billed endpoint allocation**. Peak endpoint
RSS is unavailable in these receipts; batch peaks are reported separately.

| Mode | Core / metal | MPI ranks | Wall s | Rank-wall s | Population analysis s | SCF cycles |
|---|---|---:|---:|---:|---:|---:|
| CPCM, no MBIS | qm33 La | 57 | 244.574150 | 13,940.727 | 0.385 | 11 |
| CPCM, no MBIS | qm33 Ca | 57 | 220.321270 | 12,558.312 | 0.371 | 10 |
| CPCM, no MBIS | qm36 La | 112 | 175.566832 | 19,663.485 | 1.740 | 11 |
| CPCM, no MBIS | qm36 Ca | 112 | 164.898242 | 18,468.603 | 1.213 | 12 |
| CPCM + MBIS | qm33 La | 16 | 695.604882 | 11,129.678 | 542.740 | 11 |
| CPCM + MBIS | qm33 Ca | 16 | 681.738333 | 10,907.813 | 550.388 | 10 |
| CPCM + MBIS | qm36 La | 16 | 792.605190 | 12,681.683 | 611.576 | 11 |
| CPCM + MBIS | qm36 Ca | 16 | 792.545796 | 12,680.733 | 618.369 | 12 |
| Vacuum + MBIS | qm33 La | 16 | 882.254154 | 14,116.066 | 574.688 | 37 |
| Vacuum + MBIS | qm33 Ca | 16 | 791.740796 | 12,667.853 | 587.284 | 24 |
| Vacuum + MBIS | qm36 La | 16 | 930.005890 | 14,880.094 | 707.185 | 21 |
| Vacuum + MBIS | qm36 Ca | 16 | 914.972951 | 14,639.567 | 747.686 | 16 |

Population analysis accounts for 81.6–85.1% of the internal SCF module in the
old CPCM+MBIS runs and 68.0–84.6% in the new vacuum runs. The vacuum endpoints
also needed more electronic SCF cycles; no single cause is assigned to the
whole elapsed-time difference.

Summed rank-wall seconds per pair are:

| Mode | qm33 | qm36 |
|---|---:|---:|
| CPCM, no MBIS | 26,499.039 | 38,132.088 |
| CPCM + MBIS | 22,037.491 | 25,362.416 |
| Vacuum + MBIS | 26,783.919 | 29,519.661 |

The longer MBIS wall times therefore do not translate directly into larger
rank-wall totals than the high-rank historical runs. These are observed
execution arrangements, not evidence that MBIS makes the method cheaper.

## Actual Slurm allocations and memory

Retrieved completed accounting records with `sacct`; raw output is pinned in
the audit artifact. CPUTimeRAW equals allocated CPUs × elapsed seconds, not
measured utilization. Batch MaxRSS is the accounting peak for that batch step,
not an endpoint peak or a whole-node memory requirement.

| Job | Recorded workload | Node | AllocCPUS | Wall s | Allocated core-s | Batch MaxRSS KiB |
|---|---|---|---:|---:|---:|---:|
| 1197144 | qm33 pair plus four C5B120 endpoints | node-344-8t-1 | 344 | 345 | 118,680 | 43,129,868 |
| 1196819 | qm36 La/Ca plus apo/water | node-112-1500g-1 | 112 | 513 | 57,456 | 19,738,084 |
| 1198934 | Six CPCM+MBIS endpoints, including two canonical fixed-core endpoints | node-64-768g-5 | 64 | 1,491 | 95,424 | 10,318,208 |
| 1199949 | Four vacuum+MBIS partition endpoints | node-64-768g-10 | 64 | 932 | 59,648 | 8,010,252 |

No baseline allocation is an isolated, identically scheduled four-endpoint
counterpart to job 1199949. The two 64-CPU jobs used different hosts and
different task/concurrency schedules; CPU model details are not in the
endpoint receipts. Ratios of these batch totals would mix workloads.

## Traceable artifacts and implementation hashes

The machine-readable audit copies the 12 complete execution receipts, pins
their outputs/runtime inputs/XYZ files, retains unrounded timings, and stores
the raw `sacct` result and command:

- `workspaces/global_electrostatic_20260916/cost_reference_audit_v1/audit.json`
  SHA256 `078df0c8ff653c1ff4de092239d931723f0da6b5c340e9fa1c66e5cd245dc7bd`.
- No-MBIS qm33 receipts:
  `diagnostics/pqq_boundary_pair_20260914/mxaf_1H4I_qm33/sp_mxaf_1H4I_qm33_{La,Ca}.out.execution.json`.
- No-MBIS qm36 receipts:
  `workspaces/mxaf_1H4I_qm36_qm/sp_mxaf_1H4I_qm36_{La,Ca}.out.execution.json`.
- Old MBIS receipts:
  `workspaces/affordable_challenger_20260915/pilot/1h4i_{qm33,qm36}/{La,Ca}/endpoint.out.execution.json`.
- New MBIS receipts:
  `workspaces/global_electrostatic_20260916/partition_tasks_v1/1h4i_{qm33,qm36}/{La,Ca}/endpoint.out.execution.json`.

Recorded code hashes (historical execution identities, not a requirement to
replace current files):

| Artifact | SHA256 |
|---|---|
| ORCA 6.1.1, all runs | `38b5f057452fef275c0a1b98d270ad03d0411dab70453820d1c678bd732f6c83` |
| Common task runner, jobs 1197144/1198934/1199949 | `b9278e74e5317e859e8cb05541e1ec57cb335b26dcdf65a83b89675407343b9a` |
| Common runtime renderer, same jobs | `4d75904eca4296abd892958df700063e272f5f029a0ce256bcf1b579768aeadb` |
| Queue task runner, qm36 job 1196819 | `39fad15a8f71dfa89eb56ce07c831b5f93df984a75ee4abe522502fb18a3f212` |
| Queue runtime renderer, same job | `bb3bad13b22bd8ea94a0f6350cf09acf3ebfca0d31d313266600f80b5d4ddff2` |

## Interpretation

The archived baseline provides a real cost scale, and requested population
analysis is a substantial wall-time component of the MBIS runs. A controlled
overhead factor remains **unavailable**. Total proposed production cost must
also include full-environment preparation, ESP checking and solvent solves;
their cost is not filled with zero or inferred from the quantum receipts.
This task performed receipt parsing and accounting arithmetic only. No new
calculation, benchmark variant, model or job was run.
