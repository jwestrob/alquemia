# Fixed half-step: joint improvement survives native evaluation

The single declared candidate completed all 56 native programs across 46 fixtures.
All 25 original energy/force gates pass, both paired component energy errors fall,
and all four original objective blocks improve. The small direction predicted by
the saved parameter derivatives is useful in the full nonlinear model.

This does not qualify the interaction model. Component errors remain large, and
the consumed trimer regression has not been rerun at this point. The original
SLSQP run remains failed; this is a separately declared development experiment
using an independently verified feasible half step, not an optimizer-convergence
claim or a fresh validation test.

## Component errors

P is the sum of literal metal-water interactions. B=C-P is the remaining full-shell
interaction, including all nonpair orders. These are paired Dy-minus-La errors
at fixed eight/nine-water geometries, in the admitted common-basis convention.
They do not uniquely separate polarization or three-body energy. Absolute-parent
total-energy reference components remain unavailable.

| RMS error, kcal/mol | Origin | Linear prediction | Native result |
| --- | ---: | ---: | ---: |
| Pair P | 26.682413 | 25.588723 | 25.588991 |
| Collective remainder B | 22.303152 | 21.388944 | 21.388806 |

Both RMS errors decrease by about 4.1%. Maximum disagreement between the four
predicted and actual signed component errors is 0.000283 kcal/mol.

| Signed paired error, kcal/mol | Origin | Native |
| --- | ---: | ---: |
| P8 | 25.054235 | 24.026907 |
| B8 | -20.498014 | -19.626730 |
| P9 | 28.216797 | 27.061053 |
| B9 | -23.972747 | -23.016375 |

The old E/F, pair-moment, full-shell moment and susceptibility squared norms all
decrease under exact comparisons. Their raw scales differ; they are engineering
blocks, not independent error bars. Pair-moment improvement is especially small
(5.73e-7 in its squared norm); it is not a meaningful dipole repair by itself.
No acceptance threshold was changed after this result.

## Execution and interpretation

The unchanged qualified native binaries, scientific evaluator and targets were
used on the explicitly authorized login host with twelve single-thread workers.
Wall time was 12.217 seconds; controller CPU 12.417 seconds; child user/system CPU
11.219/13.917 seconds, total measured CPU 37.553 seconds. No scheduler allocation,
new reference calls, geometry generation or held-out outcomes were involved.
Static terms and the admitted origin were reused. Zero-call preparation and
21 narrow wrapper controls passed; source review bc9116ee passed. Independent
actual-input review706c7716 verifies all46cells/56receipts/eightexports and
independently reproduces all25 physical gates. Comparison reviewc096660f
replays literal components and all four old blocks; maximum arithmetic
difference is3.55e-15. Both pass. Their CPU costs are5.188 and0.728seconds.

The finite local feasibility test is complete. The next substantive experiment
is one bounded component-aware development fit, with explicit physical acceptance
rules and a meaningful recovery criterion. The prior driver reported the 25 gates
but selected by its incomplete objective alone. Repeating tiny continuations or
changing a MACE checkpoint would not address that omission. See STRATEGY.md and
FIT_IMPLEMENTATION_OUTLINE.md for the prospective decision tree and source map.

The existing baseline, failed solve, old fit, and prior failed transfer results
are preserved. Fresh validation, sampling and production remain unreleased.
