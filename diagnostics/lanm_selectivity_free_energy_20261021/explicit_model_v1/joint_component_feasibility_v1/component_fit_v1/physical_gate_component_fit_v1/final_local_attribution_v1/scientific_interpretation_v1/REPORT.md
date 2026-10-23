# Final-point local restriction: test the nonlinear contrast first

**Recommendation:** execute only the separately declared B-half/D-half diagnostic after its source/release gates. It can distinguish a real local mechanics tradeoff from an inaccurate affine prediction. These saved results do not justify wider parameter bounds, relaxed physical gates, another fit, or a new functional form yet.

## What the saved evidence supports

At the actual final model, P/B RMS errors are **12.7203/10.4418 kcal/mol**, q=0.22727182, with all 25 original physical checks passing. The exact-final derivative panel is mechanically admitted separately by review `44bf68c0`: 24 vectors/672 programs. The following are **affine full-step predictions**, not measured candidates; the field named `actual_component_q` in the affine result is recomputed from that affine model.

| Fixed case | P/B RMS, kcal/mol | q | Scope of prediction |
|---|---:|---:|---|
| A: original box and 25 gates | 11.7433 / 9.8159 | 0.1936990 | Electrical losses reported |
| B: A plus three electrical norm protections | 11.3602 / 9.9444 | 0.1988024 | Pair/C/chi norms nonincreasing |
| C: remove only box faces | 11.7206 / 9.7970 | 0.1929535 | Outside original box |
| D: remove only physical gates | 0.01446 / 0.00313 | 0.000000294 | Fails 17 of 25 physical checks |

The strongest tested local restriction is the original energy/force block. D's paired water-minus body-force RMS reaches **15.4888 versus limit 2 kcal/mol/Angstrom**. Removing box faces buys only about 0.000746 in q relative to A; preserving the three electrical norms costs about 0.005103. Dimer radial mechanics is a major conflict, alongside active shell-force constraints. This is not a unique attribution to one radial term or microscopic process.

A/B/C retain raw `AlmostSolved` status. C's primal residual is 2.36e-7 despite positive retained scalar witness margins; D is `Solved`. These are neither optimality certificates nor a global model-capacity test. The clipped Dy_v stencil, finite trust region and first-order approximation also limit the inference. Feasible improvement remains predicted, so local impossibility is not established.

## Why apparent total-energy recovery is insufficient

Final paired P/B errors have opposing signs: P8/P9 = +11.9361/+13.4589 and B8/B9 = -8.9241/-11.7654 kcal/mol. Their cancellation leaves smaller C errors. Final pair/C dipole RMS errors are 0.639432/0.085846 D and four-state chi RMS is 0.680444 bohr^3. The W9 C dipole Dy-minus-La sign is still wrong: **+0.127900 D model versus -0.104618 D reference**. Protecting aggregate electrical norms in B does not guarantee repair of this signed contrast.

All eight consumed trimer dipole-error norms improve over the original baseline, while seven of eight E3 errors worsen. All four paired E3 errors improve over the previous fits but remain worse than baseline, at -3.5585/-2.5520/-4.1155/-2.4133 kcal/mol. Energy and electrical improvement therefore cannot be inferred from each other. Conditional chi increments have much smaller errors than the four individual states; they do not qualify the full response tensor.

The current family already has separate pair/collective drive and stiffness, plus four metal-O radius/epsilon parameters. Repeating a generic response-separation proposal or an R/epsilon-only compensation test would miss that history. A field-independent pair term cannot change fixed-geometry dipoles/chi and cancels from B and E3. A many-body energy can change collective energies but still cannot repair fixed-geometry electrical observables by itself. The current evidence selects neither modification.

## Smallest decisive next step

Root plan `06478f92` freezes **two points**, at most 112 native programs, no optimizer/reference/heldout calls and no center/probe repeats:

1. **B-half:** predicted q=0.20887175, P/B=12.0403/10.1931, all 25 gates and three electrical norms preserved. Check these actual quantities and the signed W9 contrast without automatic fallback.
2. **D-half:** predicted q=0.05681818, P/B=6.3602/5.2212, with 13 physical failures. This is a diagnostic point, never an accepted candidate merely because q falls.

If the measured contrast survives, it supports a local tradeoff between component improvement and mechanics in this family; B demonstrates only modest feasible progress. If either direction departs materially from its prediction, retain nonlinear/secant limitations before drawing a model-form conclusion. Unexpectedly feasible D would weaken the claimed obstruction. Native failure supplies no missing physical metrics. Report every gate and electrical error; no tolerance relaxation or follow-on fit is implied. No further model-form proposal is warranted before these two results.

## Reference and scope limits

P/B use paired common-parent electronic-plus-recorded-D4 cycles; absolute D4-complete parent targets and fragment-force decompositions are unavailable. The original energy/force targets are separate ECP-PBE0-D4 cycles, and trimer E3 is electronic-only consumed regression, not holdout or an exhaustive expansion of B. Raw cluster CC tightening changes one contrast by 0.063929 kcal/mol and passes its 0.25 cutoff-sensitivity gate; it validates neither CP/P/B, forces nor electrical response. SCF, basis, frozen-f, state/SOC and method limitations remain, with no quantitative physical uncertainty bound.

This interpretation reads saved results/reports only. It does not repeat numerical admission, native calls or solver work, and changes no frozen source, model, gate or shared checkpoint. Exact authorities are in `INPUTS.json`.
