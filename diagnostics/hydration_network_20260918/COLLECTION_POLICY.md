# Optimization output collection correction

Inspection of the completed native gas-water Opt output establishes that ORCA
6.1.1 performs analytic gradient evaluations during optimization, then a final
energy-only evaluation at updated coordinates. No final `.engrad` is emitted.
The final gas energy and preceding gradient energy differ by 1.81e-10 Eh and
belong to slightly different coordinates; they are not silently paired.

The original eight protein tasks mistakenly asked the shared runner to expect
an `.engrad` artifact. Preserve their manifests and running calculations. The
runner may mark an otherwise converged, normally terminated optimization failed
because that artifact is missing. Collect it with the corrected live
`hydration_network.py` into a new collection file. This explicitly retains the
receipt's missing-artifact status while separately checking scientific completion.
No calculation is rerun just to produce an unnecessary file.

The corrected collector checks normal termination, SCF and geometry convergence,
final XYZ energy/ordering, frozen coordinates (2e-5 A) and rigid water internal
geometry (2e-4 A). Printed total analytic gradients retain their own preceding
energy and coordinates, with 1e-6 A / 1e-9 Eh/bohr print resolutions. Final XYZ
must agree with the final printed coordinates within 2e-6 A. No gradient is
assigned to the later final energy-only endpoint. Unsupported numerical
derivatives remain rejected. All tolerances were declared before any protein
optimization converged.

Future prepared manifests omit the nonexistent final-gradient requirement; this
changes collection expectations only, not the quantum input. The unsubmitted
`occupancy_prepared_v1` is superseded: regenerate in a fresh directory before
execution. Do not edit the original prepared artifact or submit it unchanged.

Nine real-fixture tests pass, including the actually completed gas optimization,
reference/sign algebra, physical preparation and explicit malformed copies with
missing gradients or seeds. No fabricated scientific output is used. Missing
starts and incomplete fixed-count occupancy sets cannot become valid minima.

ORCA distinguishes Opt from the energy/gradient run type in its
[basic settings](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basics.html).
The actual completed gas output is the executable-version evidence for the
collection behavior described here.
