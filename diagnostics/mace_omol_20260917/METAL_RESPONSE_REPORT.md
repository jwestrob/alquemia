# Full3D metal response: GGR physical checks pass; alpha remains unavailable

2026-09-18. The MACE/GB curvature model predicts small GGR energy changes and
near-stationary metal positions accurately in two core representations. All four
native DFT endpoint checks pass; their score corrections differ by0.175977kcal/mol.
Alpha-lactalbumin's predicted minima exceed the frozen0.20A displacement domain,
so this protocol provides no alpha correction or demonstrated discrimination gain.
Production baseline/default and every previous experiment remain unchanged.

## Actual results

| Representation | Ca movement,A | La movement,A | Actual Ca response,kcal/mol | Actual La response,kcal/mol | Actual Ca−La correction |
|---|---:|---:|---:|---:|---:|
| GGR extended | 0.057355 | 0.081585 | −0.353519 | −1.543378 | +1.189859 |
| GGR connected | 0.048823 | 0.084584 | −0.257748 | −1.623583 | +1.365836 |
| Alpha1F6S | 0.284000 | 0.244368 | unavailable | unavailable | unavailable |
| Alpha6IP9 | 0.349366 | 0.284194 | unavailable | unavailable | unavailable |

Alpha movements are unsupported unconstrained predictions, not executed geometries.
All eight coarse/fine combined3x3 matrices are positive and pass frozen numerical
refinement checks. No negative mode clipping, trust-radius expansion or score
threshold change occurred. Actual GGR DFT+J energy errors versus prediction are
0.00118–0.03797kcal/mol; residual metal-gradient norms0.224–1.472kcal/mol/A pass
the independent gradient checks. All actual energy changes lower the center.
The GGR partition difference0.175977kcal/mol passes the2kcal requirement.

This validates a conditional small metal-response construction, with the rest of
the protein fixed. It does not validate a whole-protein relaxation, native DFT
Hessian, global minimum, entropy term or binding free energy. Long-range scaffold
electrostatics remain outside this model. Original source hydrogens are retained,
including the previously documented stretched-H limitation. No conductor correction
is combined here, and no old aquo reference or classification band is inherited.

## Why this was useful and what remains

The original single metal axis captured only1.5–31% of alpha's actual DFT metal
gradient norm. The full physical translation subspace has stable curvature and
reveals a substantial proposed displacement, while GGR's smaller response is
now independently supported. This identifies an actual limitation of the old
coordinate choice. It does not establish alpha's correct score or broad La/Ca
accuracy; both alpha−GGR comparisons remain unavailable under this protocol.

Next test: [bounded-response plan](BOUNDED_METAL_RESPONSE_PLAN.md). Apply the same
0.20A sphere to every state and evaluate a fixed finite-response descriptor.
Existing GGR solutions are interior and reusable; four additional alpha DFT
endpoints and eight short evaluations suffice. No new grid, ligand/water change,
threshold fit or radius adjustment. The old experiment remains complete and
unchanged. The next protocol has not yet been implemented or executed.

## Execution, cost and tests

Actual calls:288core MACE,224whole/core short-component evaluations,288nativeGB,
and4native r2SCAN-3c/CPCM analytic-gradient endpoints. Center results are reused
from the existing mechanics pilot; none was fabricated or recomputed here.
Jobs1201352–1201358,1201360–1201361,1201363–1201364,1201370–1201371 all completed.
Zero scientific execution failures. No new PQQ inference or baseline rescore.

Measured allocation totals: **5730GPU-seconds,127328allocatedcore-seconds,
31251.232reportedCPU-seconds**. Native DFT validation used557wall-seconds on64CPU;
matching short validation113wall/GPU-seconds. These are development costs,
not a demonstrated routine production cost. Preparation and read-only analysis
were not separately timed and are not counted as zero. Reused historical costs
remain in their original receipts. Cost files:grid_cost_v1/result.json and
validation_cost_v1/result.json in the response workspace.

Four real preparation/gradient tests pass17.195s; five additional real-grid,
fixed-native-input, actual-result, partial-qualification and corrupted-status
checks pass0.880s. No skips or fabricated successful scientific fixtures.
The GPU environment independently passed frozen minimum preflight. Initial
minimum_v1 was unexecuted because exact float equality rejected cross-environment
NumPy roundoff(max1.14e-13); minimum_v2 retains the identical predictions and
physical thresholds with a1e-10 metadata comparison tolerance.

A bounded review caught one reporting issue: endpoint validation alone must not
release a fully qualified score before the GGR partition check. This is fixed
and tested with the actual partial result. Partial outputs remain explicit,
without a baseline substitution. The original partial/final reports and their
exact reporter sources remain preserved; V2 pins its reporter inside its output.

Protocol:`DFT_anchored_MACE_GB_full_metal_translation_v1`.
Workspace:`workspaces/mace_metal_response_20260918/`.
Primary outputs:`assessment_v1/result.json`,`minimum_v2/preparation.json`,
`minimum_report_v2/result.json`. Exact commands:[runbook](METAL_RESPONSE_COMMANDS.md).

**Recommendation:** retain the production baseline and pursue the separate,
small bounded-response test. The mechanics is credible on GGR; predictive benefit
and routine affordability are still unestablished.
