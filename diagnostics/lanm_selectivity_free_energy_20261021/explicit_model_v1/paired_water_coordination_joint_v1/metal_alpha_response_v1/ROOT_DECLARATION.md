# Consumed water-cluster response to metal polarizability

## Question and prior evidence

The terminal joint six-parameter fit passes22/25 scientific controls. Its W8→W9
paired energy error remains−4.338284kcal/mol; all parameters are interior and the
last full trial batch produces no improvement. The selected-point linear response
has almost no remaining least-squares descent. This establishes a local result,
not global infeasibility. A separate saved linear-feasibility check asks whether
the unchanged individual gates could pass at a point with higher aggregate loss.

This fixed probe asks whether the previously fixed metal polarizabilities provide
additional parameter-response directions relevant to the residual coordination
error. It does not fit polarizabilities or establish their physically correct values.
The current AMOEBA model already contains mutual polarization. In its native
implementation changing alpha changes both susceptibility and the damping length;
these effects are not uniquely separated by this experiment.

## Frozen inputs and finite matrix

Use the exact terminal17676 selected geometry, six fitted parameter values, native
engines and ten consumed reference cells, admitted by selected-ten review17683.
The expected scientific flag is false; exact raw/score review must pass. Do not
use the separate fresh-validation gate or expose fresh radial/orientation results.

Metal polarizabilities remain anchored to La1.6597 and Dy0.7280. Perturb one metal
at a time by relative−5%,−1%,+1%,+5%. For each perturbation evaluate its five
consumed fixtures: water_minus,water_central,water_plus,aquo8_central,aquo9_central.
This is8 probe vectors,40 new cells and80 native programs (analyze plus analytic
Cartesian testgrad for each cell). Reuse the other metal's exact five selected
cells when assembling each full ten-cell diagnostic. Keep all geometries, states,
reference energies/forces, charge+3, multipoles, other polarizabilities, Rmin,
epsilon, damping numbers and all other settings unchanged. Change only the numeric
alpha token on the matching metal POLARIZE line. Record exact token and resolved
parameter checks before execution. Preserve original source key bytes.

One standard_expand task,8CPUs/default memory/time; at most8 single-thread native
workers, hard ceiling80 programs. No Cartesian finite differences, quantum calls,
new geometry, optimizer, automatic retry or adaptive extra probe. Global stop on
native/admission failure; retain every missing/failed slot. Perform source/prepared
input controls and independent review before root release.

## Fixed analysis

For each eight full probe matrices report the unchanged413-residual objective and
all10 energy/15 force checks, retaining the baseline and both signs/scales. These
are diagnostic responses, not candidate selection. Compare central relative-alpha
finite differences at1% and5%; report scale dependence explicitly.

Join both extra derivative columns to the saved six-column selected-point response
matrix using the same residual units and1e-6 relative singular cutoff. Report rank,
remaining residual, predicted W9 energy error and implied parameter changes at
both scales. Scaling the two new coordinates means relative alpha change; retain
this convention alongside the original six normalized coordinates. A projection
is unevaluated and unconstrained unless explicitly checked; it cannot qualify a
candidate. Distinguish nonlinear-scale response from numerical failure.

Preserve alpha/damping coupling, provisional ECP-reference limitations and the
prior17670 review-loss notice. Additional linear response would motivate a separate
physical/parameter qualification; absent response would limit this local direction.
Neither outcome proves a unique missing term or global impossibility. No fit,
combined donor/water model, fresh validation or bulk/protein release follows.
