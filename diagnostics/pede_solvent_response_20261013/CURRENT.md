# PedE solvent-response checkpoint — 14 October 2026 UTC

The subsequent [inconclusive-source review](../pqq_inconclusive_review_20261014/REPORT.md)
is complete:48 peers/nine ambiguous contexts, no new chemistry or score changes.
Read that checkpoint for the current diagnosis and proposed physical next test.

**Latest: the fixed-charge candidate completed and corrects both PedE errors.**
All five two-Ca PedE sources are now correct. The recovered125 source set changes
119correct/3wrong/3inconclusive to121/0/4. However, five previously correct
three-source summaries become inconclusive; no blanket promotion. Read
[the candidate report](frozen_distribution_candidate/REPORT.md) and EVIDENCE.json.
All candidate jobs4006–4010/4014/4016 are terminal; do not repeat them from delayed wakes.
Production and all historical results remain unchanged.

**All eight DFT endpoints and comparison3377 are complete.** The source-dependent
Ca−La solvent contrast is +10.0342163383 kcal/mol with r2SCAN-3c/CPCM and
+21.2566283014 with corrected GFN2/ALPB. The direction survives; the magnitude
and individual metal responses are strongly model-dependent. The DFT comparison
alone did not establish a classifier fix; the later candidate test above is separate.

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

## Native fixed-charge replay and candidate complete

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

The [research candidate](frozen_distribution_candidate/PLAN.md): MACE vacuum
energy plus native ALPB evaluated on its matching GFN2vacuum charges. Keep the
full continuum functional, explicitly omit extra solute electronic relaxation,
and test usefulness with separate canonical calibration followed by frozen-band
structural transfer. No fitted response fraction or production promotion.
The exact inventory and calibration rule were frozen; canonical4007/4008 passed
before transfer4009/4010. All1542 native inputs are available, with3076 new component
calls and four reused input pairs. No new MACE/SCC/DFT or optimizer work occurred.
Canonical25 remains25correct; expanded135 changes134/1/0/0 to129/0/6/0. Existing
100triples change94correct/6unavailable to89correct/5inconclusive/6unavailable.
The genuine PedE gain and abstention regressions remain separate. Candidate test,
replay/calibration, summary and saved-matrix diagnosis cost731 allocatedCPU-s,
zeroGPU, excluding previously computed MACE energies/charges. Full phase total
1,077,214allocatedCPU-s.

Saved-matrix diagnosis4016 also completed: new abstentions arise from energy-model
changes plus stricter Ca bands, not harmful geometry reselection. PedE sample2's
new correct margin is only0.233kcal/mol. Read selection_diagnosis/REPORT.md inside
the candidate directory. Next is a physical response-model decision; no next
molecular run or promotion is queued. Do not refold or retune bands to rescue it.
All scientific computation must continue through Slurm with completion wakes.

## Queue and recovery

All PedE DFT/comparison/replay/candidate jobs are terminal. Root read/acknowledged
3372/3377,4000–4002,4006–4010,4014 and4016. Delayed wakes need acknowledgement, not reruns.
Scheduler still lists unrelated reporting2669/2832/2833 on node-48-256g-18 despite
completed output receipts; preserve these under Jacob's instruction. Readiness
probes remain off. Candidate agents completed inventory/implementation/independent
review; root owns execution, shared scoring/status and final delivery.
