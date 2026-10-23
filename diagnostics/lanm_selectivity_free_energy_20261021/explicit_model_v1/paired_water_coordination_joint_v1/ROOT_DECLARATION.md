# Joint water coordination development test

Declared 21 October 2026, after the completed frozen-candidate W9 screen and
before any new parameter evaluation. Root owns integration and submission.

## Question and evidence

Can the existing six-parameter metal–water family meet the same energy and rigid
force limits simultaneously for the consumed dimer radial fixtures and both
eight- and nine-water configurations?

Job17670 failed three of four coordination/geometric transfer checks. Its paired
W8→W9 energy error is −4.951230 kcal/mol. The saved original model error is
−5.893211 kcal/mol: the water fit improved the energy contrast but worsened the
individual W9 forces while improving W8. These facts motivate a joint development
objective. They do not establish a unique missing physical term or prove that
the parameter family cannot represent both configurations.

## Frozen inputs and allowed changes

Use the exact eight cells and consistent provisional uncorrected ECP-PBE0-D4
surface from the17605 water fit, plus the two original W9 cells admitted in the
17670 preparation. Preserve every geometry, charge, state, reference output and
failure. W9 also rearranges the original eight waters; do not call this a pure
ninth-water insertion. No new quantum calculation is included.

Parameter order remains La radius, La damping, Dy radius, Dy damping, La epsilon,
Dy epsilon. Change only metal–water-O VDWPR radius/epsilon and POLPAIR damping.
Use the exact original six-parameter intervals and native engines; do not widen
bounds or change multipoles, polarizabilities, water-H, donor or bonded terms.
Start from the exact frozen17605 selected vector. Its eight old model cells and
both17670 W9 cells must replay at zero new calls before execution.

All ten configurations are now development data. The four previously exposed
orientation results remain historical validation/regression evidence, not a new
blind test. They are excluded from fitting. No reserved donor or protein outcomes
may be opened by this preparation or run.

## Joint objective and unchanged individual acceptance limits

Keep four equally weighted blocks: energy, La rigid response, Dy rigid response,
and Dy-minus-La rigid response. Each block has weight0.25.

The energy block gives weight0.5 to four aquo8-referenced paired contrasts:
the three dimer radial configurations plus W9, with equal weights within that
subblock. Give the other0.5 to the four existing individual radial-work residuals,
again equally weighted. Retain the two paired radial works as acceptance checks.
This gives eight weighted energy residuals and ten separate energy controls.
Every energy control must satisfy absolute error≤1 kcal/mol.

Each force block averages the five geometry fixtures equally, using the existing
category-balanced net-force and centroid-torque/1Å projection unchanged.
Every one of the15 individual fixture/category checks must satisfy rigid RMS≤2
and maximum component≤5 kcal/mol/Å. No average can hide a failed fixture.
The joint residual vector has413 entries in the declared geometry/category order.
Report old eight-cell controls, added W9 controls and all residual blocks separately.

## Finite algorithm and resources

Reuse the existing six-dimensional algorithm: clipped normalized±0.005 probes,
least-squares cutoff1e−6, componentwise normalized trust box±0.025, scaled
Gauss–Newton and bounded least-squares proposals at half/full length in the same
fixed order. Complete every derivative/trial batch before selection; current
point wins ties and derivative probes never become the selected point.

One eight-CPU standard_expand task, at most eight single-thread native workers,
default memory/time. Maximum24 outer steps,384 new full vectors (385 including
the initial vector), and5120 new native program invocations. A full step needs
at most120 derivative and80 trial programs before exact effective-input reuse;
24 steps therefore fit within the declared ceiling. No new Cartesian FD is
included. Molecular execution follows actual preparation, meaningful saved
controls and independent source/input review.

Stop on all individual controls passing, no actual improvement after a complete
trial batch, an explicit ceiling, or native/admission failure. Preserve the actual
stop, parameter-bound contacts, residual span and counts. No automatic extension,
changed weights, tolerance, parameter family or second starting point is included.

## Interpretation and coupling to other work

A pass shows representability on this expanded consumed development set against
the provisional reference. It does not validate physical hydration, salt,
periodic dynamics, protein transfer or exchange free energies. Register fresh
validation before claiming renewed independent transfer evidence.

A local stop or bound contact is diagnostic, not proof of global impossibility.
Keep the original water candidate and all failed tests immutable. Donor fit17621
continues with its original water model and settings; its result cannot silently
be combined with changed water parameters. Any later combined model requires
matched wet/dry donor reevaluation and the existing numerical/physical gates.
