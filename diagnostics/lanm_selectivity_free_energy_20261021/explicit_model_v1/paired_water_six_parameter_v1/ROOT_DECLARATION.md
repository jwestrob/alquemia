# Bounded six-parameter water development

## Question and decision

Can joint variation of metal–water pair radius, polarization damping and well
depth reproduce the consumed water energy and force targets under the existing
criteria? The four-parameter fit17593 reached a local stop. Diagnostic17596 adds
two independent well-depth directions, with parameter-linearity checks passing.
Its nested linear projection changes weighted loss0.4194902375 to0.1403437454.
This is motivation for a new actual fit, not an achieved improvement or a proof
that a physically adequate model exists.

## Frozen data and family

Use the same eight consumed water cells, training references, weighted objective,
per-fixture RMS≤2 and maximum-component≤5 kcal/mol/Å, and energy-error≤1 kcal/mol.
Keep all historical outcomes. The initial four parameters are the exact selected
17593 vector [3.832755583060068,0.31164088100092946,
3.525530001514476,0.32135326991694113]. Their original radius/damping bounds stay.

Add only La and Dy explicit metal–water-O VDWPR epsilon, initially
0.24429636467571708 and0.37332330078584813 kcal/mol. Freeze each new interval at
0.5–1.5 times its original value. These positive bounds define a development
experiment; they are not estimated physical uncertainty. Other pair terms,
water-H interactions, multipoles, polarizabilities and bonded terms stay fixed.
No geometry, electronic reference, solvent correction or target is changed.

## Algorithm and finite cost

Retain the existing finite-difference/bounded-trial algorithm: each parameter is
normalized to its frozen interval; central probes use±0.005 with recorded clipping
and actual widths; trust radius0.025; the same four half/full Gauss–Newton and
bounded-least-squares trial proposals. Select only actual objective decreases;
current point wins ties. At most24 complete steps,384 new vectors and4096 native
program invocations. Twelve derivative vectors plus four trials per complete
six-parameter step give the declared384-vector ceiling. Reuse exact admitted
physical cells and avoid recalculating an unchanged metal. Preserve full cell
availability and native provenance. No automatic retry or new numerical settings.

Stop on passing all original controls, a complete trial batch with no actual
improvement, the finite budget, or native/admission failure. Bounds reached must
remain visible. Do not expand bounds, retune thresholds or switch objectives
within this experiment. Run on one eight-CPU standard_expand allocation, up to
eight single-thread workers, default memory/time; root owns submission and wake.

## Interpretation and subsequent work

Inspect actual energies, forces, residual blocks, displacement of parameters,
stopping reason and costs. Passing these consumed fixtures permits separately
declared transfer tests; it is not independent physical validation or an exchange
free energy. No donor fit, reserved validation, protein simulation or production
promotion is part of this run. Preserve diagnostic17596 without selecting a probe
as a fitted candidate. Preparation and meaningful saved controls precede release.
