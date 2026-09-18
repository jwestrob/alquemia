# Completed: unchanged typed-group PQQ model

**Final status:** job 1201524 completed all 50 endpoints. All 54 charge checks
pass; calibration fails with class overlap and 106/154 expected orderings.
See REPORT.md, RESULT.json and VALIDATION.md for final results and costs.
The following initial-execution snapshot is retained as history.

## Historical running snapshot

Job **1201524** is running on one A5000, 16 CPUs and 64474 MiB host RAM.
Manifest V2 SHA:
`3509b45a52618c0e61a79ceb7c76db180f06a4f63f36efa1b0d40a5575705192`.
Do not resubmit it. Scientific protocol and production default are unchanged.

At `execution_progress_v1.json`, six of 50 new endpoints had accepted receipts,
with no failed molecular evaluations. Inference took 61.38–63.73 seconds per
endpoint; summed inference 377.402845 s, peak GPU allocation 10,568,285,696 bytes.
These are partial measurements, not final job cost or calibration results.
The first three pairs comprise two La references and one Ca reference; all
remaining cases must finish before applying the frozen calibration rule.

Implemented: `scripts/mace_group_canonical.py` uses the existing typed group
builder, endpoint construction, charge adapter, runner, locks and receipt cache.
It pins original labels and all 25 full preparations, reuses four identical
crystal endpoints, retains missing 1KB0 and prior grouping failures, and reports
uncalibrated contrasts alongside the historical baseline. No new threshold is
available while the run is incomplete. No new solvent, DFT or training calls.

Five real source/parser/calibration/cache tests passed in 238.311 s, zero skips.
These use actual baseline and earlier failed OMOL scores as algebra fixtures,
never as new MACE evidence. A further actual-partial-report test passed in
0.003 s, zero skips: missing endpoints remain null, all classes/bands unavailable.
New complete molecular integration and final report tests remain pending.

Group preparation took 77.917176 wall and 76.674563 CPU seconds. Model V1 failed
local preflight because two transitive source modules were absent from its
snapshot, before any molecular call/submission. V2 adds those files; exact
model, software, task settings and all 50 coordinate hashes are unchanged.
Both versions and the failure are preserved. Other local preflight/test/report
CPU time is additional and not fully profiled. Final allocation costs pending.

See [plan](PLAN.md), [source correction](SOURCE_CLARIFICATION.md) and
[commands](COMMANDS.md). The correction identifies 1KB0's actual missing-peptide
source problem, replacing the handoff's erroneous TRO explanation. Full source
pins, execution receipts and partial reports live under
`workspaces/mace_group_canonical_20260918/`.

Next: allow the fixed run to finish, collect/report all 25 calibration cases and
three transfers, record actual costs, and test the complete output. Bands are
allowed only for a complete positive gap; otherwise report calibration failure.
Neither outcome removes the GGR grouping warning or authorizes promotion.
