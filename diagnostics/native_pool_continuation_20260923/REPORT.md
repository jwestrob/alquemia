# Native continuation settles these pools without improving their calls

**All four fixed pools pass the declared numerical settling test; their selected
geometries and classifications are unchanged.** The persistent A0A3F2YLY8
Ca-conditioned sample 1 error therefore is not resolved by further native SCF
settling of these existing candidates. This closes this four-pool numerical
test; it does not justify a larger accuracy experiment by itself.

Both passes completed all 80 cells with actual XTBRESTART/native-mixer evidence,
unchanged states/parameters and valid receipts. Stage 2 supplies all eight
requested analytic gradients to the separate force-check experiment. No cells
failed or were substituted, and no third continuation was run.

## Discriminatory result

The reported endpoint is always stage 2, as frozen before execution. Larger R
is more La-like. Differences below are model kcal/mol; all unrounded components
and all 20 same-geometry contrasts remain in the pinned result.

| Consumed source | Stage 2 R | ΔR from archived pool | Adaptive old-band call, before → after |
|---|---:|---:|---|
| 1H4I | −405505.439970114 | −0.019919558 | Ca → Ca, correct |
| 4MAE | −405444.287782274 | +0.016643087 | La → La, correct |
| Q88JH5 | −405464.118096698 | +0.004599565 | Ca → Ca, correct |
| A0A3F2YLY8 Ca-sample1 | −405463.300765336 | +0.004753490 | Ca → Ca, wrong |

Coverage is identical: four complete pools, 20 geometries and 80 cells per
stage. Adaptive-band transfer remains **3 correct / 1 wrong**. Released-static
band transfer remains **2 correct / 2 inconclusive**, with Q88JH5 and A0A3
inconclusive. These are diagnostic transfers through existing bands, not a
newly calibrated classifier or independent affinity validation. Own calibration
is unavailable; none was fitted to these four cases.

All eight selected metal-row geometries remain unchanged. The mathematical and
operational pool selections agree here; the existing 0.1 kcal/mol origin
selection rule was preserved. No alternate geometry, stage or threshold was
chosen to improve a label.

## Numerical result and its limits

| Stage 1→2 check | Maximum observed magnitude | Frozen limit |
|---|---:|---:|
| Individual native cell energy | 0.012678915 | 0.1 |
| Same-geometry composite Ca−La contrast | 0.017007658 | 0.2 |
| Pooled composite Ca−La contrast | 0.004723191 | 0.2 |

All 80 cell checks, all 20 same-geometry checks and all four pools pass.
The largest stage 2 shift from any archived native cell is only 0.041296003
kcal/mol. Thus the earlier approximately 4.805 kcal/mol Q88JH5 near-identical
**new-template** anomaly is not reproduced by these archived five-candidate
pools. The earlier confirmed restart repair remains useful evidence about that
specific anomaly; it was not evidence that these existing pools were defective.

Passing two continuations establishes this declared local energy-stability
test. It does not establish a global electronic minimum, universal native SCF
accuracy, force accuracy, equilibrium populations or biological affinity.
Printed generic SCF residuals are retained separately in each row. The
80 stage 2 outputs all satisfy their printed energy-change criterion but all
exceed the printed generic MAX/RMS density criteria; these displayed criteria
are not claimed as satisfied native stopping conditions. The
derivative-check agent will evaluate the shared eight actual gradients; no force
accuracy claim is made here.

## Execution, integrity and cost

Job **1210333** completed 160 native GFN2 calls, including eight analytic-gradient
calls, in **381 seconds with 64 allocated CPUs: 24,384 allocated core-seconds**.
Allocation: 128 GiB on node-224-2t-8gpu-1, CPU-only gpu partition, zero GPU time.
Stage executor wall times were 177.111159 and 170.229965 seconds; scheduler time
also includes in-job preparation/collection. Local setup/reporting is additional
and unmetered. Slurm batch MaxRSS was 5,031,768 K; this is the scheduler's reported
step statistic, not an independently measured whole-allocation peak. No new
MACE, DFT, geometry optimization, protonation or water changes occurred.

All pre/post GBW/xtbw files, raw outputs and receipts remain available. Six
focused actual-artifact/parser tests pass, zero skips after completion; prelaunch
the scientific completion test was explicitly unrun. This testing is separate
from the 160 actually executed molecular calls.

The running implementation remained immutable. A report-only correction separates
archived `low` receipts from the continued energies in the final matrix. The
original automatic `result.json` is retained; **`result_receipts_v2.json`** is
the authoritative comparison. Every score, decision, check and selection is
identical between those versions. Only receipt attribution changed.

Baseline/default and all historical outputs are unchanged. Recommend retaining
the released scorer while the separately approved derivative/scaffold work
tests other explanations; do not widen this continuation experiment merely
because it produced no classification repair.

Exact artifacts: [ARTIFACTS.json](ARTIFACTS.json). Runnable no-compute reporting
and test commands: [COMMANDS.md](COMMANDS.md). The fixed scope is in
[PLAN.md](PLAN.md).
