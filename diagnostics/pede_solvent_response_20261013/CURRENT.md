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

## Native fixed-charge replay complete; next candidate preparing

[The native replay](fixed_charge_replay/REPORT.md) completed as4002: all eight
fixed-charge evaluations and16water-component gates pass. The original21.2566
kcal/mol source/metal contrast separates into2.5099 on vacuum charges plus
18.7467 from additional redistribution including intrinsic compensation. The
closure residual is1.1e−5kcal/mol, within the predeclared precision budget.
This identifies model sensitivity, not a proven erroneous contribution.

Two technical failures are preserved:4000CMake export,4001wrong wrapper XYZ
conversion import. Version3 uses the actual reader conversion; native parameters
and all tolerances stay unchanged. Final build/evaluations65s on12CPUs; all three
attempts1668allocatedCPU-s, zero newSCC/DFT, zeroGPU. Entire current phase through
4002 totals1,076,483allocatedCPU-s. All actual results/costs are in EXECUTION.json.

Next [research candidate](frozen_distribution_candidate/PLAN.md): MACE vacuum
energy plus native ALPB evaluated on its matching GFN2vacuum charges. Keep the
full continuum functional, explicitly omit extra solute electronic relaxation,
and test usefulness with separate canonical calibration followed by frozen-band
structural transfer. No fitted response fraction or production promotion.
Inventory/code preparation runs in parallel; no candidate evaluation submitted.
Canonical rule/population must be pinned before execution. All new compute uses
Slurm and actual completion wakes. Do not rerun completed DFT or replays.

## Queue and recovery

All PedE DFT/comparison/replay jobs are terminal. Root read/acknowledged3372,
3377result and4000/4001/4002. Delayed wakes need acknowledgement, not reruns.
Scheduler still lists unrelated reporting2669/2832/2833 on node-48-256g-18 despite
completed output receipts; preserve these under Jacob's instruction. Readiness
probes remain off. New candidate inventory/preparation agents are active; root
owns execution and shared scoring/status.
