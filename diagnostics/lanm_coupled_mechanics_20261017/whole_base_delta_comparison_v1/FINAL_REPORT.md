# Electronic-base candidate: closed after development test

19 October 2026. **The candidate failed its predeclared comparison and will not advance.**
Four 7CCO development states are complete; four 7CCN candidate states remain unrun.
The independent 7CCN PBE0 reference job14795 continues under its original gates.

The tested energy was whole GFN2 + regional MACE − matched regional GFN2 +
existing fixed-charge ALPB. All structures, states and La/Dy parameters were fixed.

## Result

Deformation cost, kcal/mol, for the same before/after coordinates:

| Model | La | Dy |
|---|---:|---:|
| Whole PBE0 + same ALPB | 12.280 | 11.723 |
| Whole MACE + same ALPB | 11.290 | 11.082 |
| Whole GFN2 + same ALPB | 10.399 | 10.794 |
| Whole GFN2 with regional correction + same ALPB | 9.871 | 9.200 |

The candidate undershoots PBE0 by **2.409 / 2.523 kcal/mol**, failing the frozen
1-kcal/mol development limit for both metals. The regional correction also worsens
the base-only work discrepancy. Whole MACE remains closer on this consumed case.

The new model reduces the old hybrid's large force errors, but is less accurate
than whole MACE here: La force-error-change RMS is2.028 versus0.530 kcal/mol/A
(old hybrid4.464). These two endpoints do not establish a Hessian, converged
minimum, transferable model or metal preference. Existing electronic-state,
terminal-chemistry, frozen-f and finite-temperature qualifications remain absent.
Missing Dy comparator forces remain missing. The local correction is not mutually
polarizable and the solvent term remains additive.

## Decision and next step

Close this candidate's four conditional 7CCN calls and optimization path. Do not
retune the gate or replace its failures. Continue the separate whole-PBE07CCN
reference to test whether whole MACE's better mechanics agreement transfers to
the second existing structure. The broader mechanics goal remains unmet.

## Reproducibility and cost

Six new molecular calls total:4wholeGFN,1regionalGFN,1regionalMACE. All other
components were reused. Saved assembly14805 and comparator14806 make zero calls.
Including the preserved pre-forward failure14803, total allocation is
**260 CPU-seconds and14 GPU-seconds**; this is not utilization. Independent DFT
is excluded from that cost.

GPU14803 failed before model load on bitwise re-creation of admitted cap-gradient
arrays. Versioned source50076ac2 reuses the pinned admitted arrays with the raw
native/input/state/map checks retained, without new numerical tolerances.
Recovery14804 passes exactly one forward. Original failure and raw outputs stay.

[Full comparison](actual_result_v1/execution/REPORT.md),
[raw result](actual_result_v1/execution/RESULT.json),
[root decision](actual_result_v1/ROOT_DECISION.json), and
[independent review](INDEPENDENT_RESULT_REVIEW.md).
Result SHA:ae68645d451db9f41dc6dda328aaae0cf94a59113aedc1c5d48d79046e919236.
