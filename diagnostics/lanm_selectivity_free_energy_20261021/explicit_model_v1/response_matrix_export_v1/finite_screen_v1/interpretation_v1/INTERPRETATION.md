# Finite drive/stiffness screen: useful response, incompatible tested tradeoffs

**Decision:** implement one conservative, regularized pair/collective response
extension together with the existing metal–O short-range energy parameters.
Do not select a grid point or repeat this grid. The screen establishes separate
drive and stiffness leverage, but its tested points do not repair the consumed
water model. The proposed implementation is specified in [SPEC.md](SPEC.md).

## What actually ran

Saved RESULT `19cb1ea24b6923c3ab5ca9c64d265d31a590b6a3ddc703e78d9f1a53d6fe0686`
contains 486 matrix states: nine points per metal on 27 fixtures each, composed
into all 81 metal-pair combinations. Independent review `dfc04b6c` passes;
all states are stable (normalized eigenvalues 0.414172–2.137701, maximum
condition number 4.375922). Thus instability does not explain these tradeoffs.
The saved screen took 3.176 s wall/2.794 CPU-s on one thread, with zero native
calls. All 15 candidate force checks remain unevaluated.

Every combination fails the original W8/W9 coordination-energy check:
errors span **−5.814256 to −2.743836 kcal/mol**, versus baseline −4.338284 and
the unchanged absolute 1 kcal/mol gate. The maximum remains 9/10 energy gates.
The counts at 0, 1, 2, 3, 7 and 9 passed gates are respectively
18, 8, 38, 8, 4 and 5 combinations. All five 9/10 combinations have zero added
drive for both metals. All 81 results, including losses, remain in FACTS.json
and the immutable producer result; no point is selected.

## Exact response tradeoffs

The table uses the already-consumed reference h=0.001 au energy estimator for
readability. Both reference steps/estimators are retained in FACTS.json. These
differences have no newly invented physical-equivalence threshold.

| State χzz (a0³) | Reference | Baseline dense | Positive stiffness, l=+0.385081767 |
| --- | ---: | ---: | ---: |
| La dimer | 15.472822 | 18.819707 | 15.455969 |
| La trimer | 22.016524 | 26.277555 | 22.337551 |
| Dy dimer | 13.158107 | 12.926079 | 11.455964 |
| Dy trimer | 19.880690 | 20.785367 | 19.070569 |

Positive metal stiffness greatly reduces La's susceptibility deviations, but
reduces its already weak pair moment. At zero drive, W9 Pz falls
2.916726→2.730086 D, against reference 3.929195 D. Negative stiffness increases
both Dy state susceptibilities, while positive stiffness decreases both. Drive
does not change χ in this linear functional. Consequently, one Dy stiffness
amplitude cannot move its dimer and trimer toward their two exact reference
values simultaneously. This is an analytic directional limitation of this
family, not an accuracy claim based on the small dimer offset.

| W9 quantity (Debye, laboratory z) | La baseline | La tested range | La reference | Dy baseline | Dy tested range | Dy reference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| C | 2.107200 | 1.846488…2.569075 | 2.244444 | 2.017553 | 1.831397…2.306344 | 2.139826 |
| P | 2.916726 | 2.563470…3.378349 | 3.929195 | 2.729250 | 2.467165…3.051058 | 3.959698 |
| B=C−P | −0.809526 | −0.872647…−0.687515 | −1.684751 | −0.711697 | −0.770236…−0.630407 | −1.819872 |

The total C can bracket its reference while pair response P and collective
suppression B remain insufficient at every tested local point. The needed
changes are approximately +1 D in P but only +0.12–0.14 D in C. A correction
that raises every water response together must therefore also gain an independent
additional-water suppression mechanism. Agreement of C alone would hide the
wrong cancellation. Symmetry-suppressed W8 total vectors do not establish accuracy.

At a=+0.1 and l=0, La/Dy Pz increases by 0.174011/0.165206 D. The same changes
lower the literal W9 water-1 dimer energies by **5.778723/6.118448 kcal/mol** and
W9 cluster energies by **28.543711/31.256428 kcal/mol**. These are endpoint
energy shifts, not independently gated observables; the original contrast/work
gates show the resulting incompatibility. The direction has response leverage,
but fixed previously fitted short-range energies cannot be assumed compatible.

Every sampled trimer E3 remains below its reference at every local point.
Across the eight contexts and all tested points the errors span
−8.983817…−2.708254 kcal/mol. All eight μ3 magnitudes also remain too small.
For example, W9 far La spans 0.046645–0.097080 D against 0.349706 D, and Dy
0.038396–0.064634 D against 0.409985 D. Independent amplitudes can change some
metal orderings by trading one error against another; that does not remove the
individual deficits. Eight E3 values are not a decomposition of the full-shell B.

## Consequence for implementation

Root's next family has **eight response amplitudes**: isolated-pair drive,
additional-water drive, pair-dependent metal stiffness and additional-water
metal stiffness, separately for La/Dy. Jointly refit only the four existing
metal–O radius/epsilon parameters; freeze POLPAIR and native damping lengths.
The extra drive direction prevents pair response from being forced to track
full-shell response; the extra stiffness direction removes the single-sign Dy
curvature restriction. These are algebraic capabilities, not proof that the
family will fit or transfer.

This differs from the failed alpha, positive EXCHPOL and D/M probes: it changes
the driving field and local inverse response independently, with environment
terms absent in an isolated pair. It preserves the original water-only
Hamiltonian. Existing EXCHPOL and D/M keywords do not implement this functional.
The new drive must be regularized at contact before native implementation;
the completed grid used an unregularized r⁻² added field and does not numerically
qualify the distinct regularized family.

Follow the frozen joint DATA_CONTRACT: original ten E/F cells and all 25 gates,
34 pair vectors, four full-shell C vectors and four independent χ states support
development. Keep all eight consumed E3/μ3 cycles outside the fit as regression;
P/B and field increments are derived reports. The separate objective declaration
must precede fitting; printing bounds are not statistical weights. Fresh W9
radial, donor-orientation and protein slots stay excluded. A finite fit may
still fail because the response family or radial energy form is inadequate,
because parameters are unidentifiable, or because provisional references conflict.

Only the tested points are rejected by this grid. It neither excludes the
continuous four-amplitude domain nor proves a unique missing physical process.
The exact Dy sign restriction is narrower and stronger. Energy deterioration
with fixed nonpolarization terms does not rule out a joint conservative refit:
those terms may compensate for the original response. Native analytic gradients,
energy-conjugate moments, stability and untouched-water controls must be qualified
before any fitting or transfer claim. No new points, fit, engine change or native
execution was performed for this interpretation.

Full paths/hashes, baseline values, all local profiles and all 81 energy results:
[FACTS.json](FACTS.json), SHA `a3889dfffb8cd57875a5c3cc1f7d3b313f89c5fd3c1b09829d19c4db8fc9130a`.
