# PedE solvent-response checkpoint — 14 October 2026 UTC

**All eight DFT endpoints and comparison3377 are complete.** The source-dependent
Ca−La solvent contrast is +10.0342163383 kcal/mol with r2SCAN-3c/CPCM and
+21.2566283014 with corrected GFN2/ALPB. The direction survives; the magnitude
and individual metal responses are strongly model-dependent. No solvent-model
fix or improved classification has yet been demonstrated.

Read [REFERENCE_REPORT.md](REFERENCE_REPORT.md), [the scientific comparison](comparison/REPORT.md)
and [electronic limits](solvent_audit/REFERENCE_ELECTRONIC_LIMITS.md).
The result does not justify removing solvent, fitting a scale factor, replacing
PQQ geometry or changing bands. The DFT reference has unresolved electronic-state
warnings and is not a ground-state oracle. Production and failed classifications
remain intact.

## Completed work and failures

Jobs3354/3371–3374 evaluated exactly eight frozen192-atom origins: two PedE
sources × Ca/La × vacuum/water. No molecular retries or optimization. The
scout's original batch failed only during parsing; recovery3370 passed15 actual
parser/fixture checks and recovered its converged result. The saved comparison
passed six partial-fixture checks in3376 and collected all eight cells in3377.
Original parsing/cofactor-reader/monitor failures remain preserved.

Final3372 completed50971s on12 CPUs. Its calcium-vacuum cell accepted TRAH
residual6.267086e−6 under the printed1e−5 threshold. Earlier descriptions of125
as a total iteration cap were incorrect: unchanged native input prints125,
then executes50DIIS iterations and TRAH through macro98, reporting843SCF cycles.
No limit override or seeded fallback was run. The prepared same-state seeded
contingency is closed; do not launch it from old commands.

Molecular cost1,074,804 allocated CPU-seconds, zero GPU; current phase including
saved audits/parser/collection jobs totals1,074,815 allocated CPU-seconds.
These are development-reference costs, not routine classifier runtime.
Exact receipts, costs and result hashes: [EXECUTION.json](EXECUTION.json).

## Next authorized diagnostic

[FIXED_CHARGE_REPLAY_PLAN.md](FIXED_CHARGE_REPLAY_PLAN.md) defines eight inexpensive
native ALPB component evaluations using the existing vacuum/water charge
vectors. It separates a frozen-distribution contribution from additional
solvent-driven redistribution while retaining intrinsic-energy compensation.
First qualify four water-charge diagonals against archived native components.
Use the actual eight-decimal charge precision and frozen analytic rounding
bounds; no fabricated shell populations, new electronic calculations or fitted
parameters. Agent la_parameter_scope owns wrapper preparation, root owns
submission/scoring. An independent agent reviews the interface and accounting.

At this checkpoint the isolated xTB6.7.1 wrapper/build is being prepared, **not
submitted**. Compilation/tests/evaluations use standard Slurm under current
shared policy, no memory/time overrides. Record new submission and completion
wake before yielding. Read the replay READY/SUBMISSION records for a newer state.
Do not rerun the eight DFT endpoints or historical CPCM/CHELPG/LanM campaigns.

## Queue and recovery

All owned PedE molecular/comparison jobs are terminal. Scheduler still lists
old reporting2669/2832/2833 on node-48-256g-18 despite completed output receipts;
preserve these unrelated entries under Jacob's instruction. Readiness probes
remain off. Completion wakes for any newly submitted work remain mandatory.
Root has read and acknowledged3372 and the actual3377result. A delayed terminal
wake requires acknowledgement, not duplicate chemistry.
