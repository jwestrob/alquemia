# Joint metal and selected-donor proposal — preparation only

Protocol: `common_metal3_four_angular_native_OMOL_SLSQP_proposals_v1`.

Parent explicitly authorized preparation on 2026-09-22; no optimization or
solvent submission is authorized by this checkpoint. Same consumed development
sources: 1H4I, 4MAE, and PLM sample0 sequences 83440678cbbd658047c9 and
07ab500e3df76b30d71c. The PLM proteins have no established biological labels.

## Motivation from the completed angular pilot

Actual final native-MACE metal-gradient norms are 13.33–40.27 kcal/mol/Å.
Seven of eight selected angular spaces have residual loads below 0.0007
kcal/mol/Å; the eighth is the already reported 0.8 Å-boundary PLM8344 La endpoint.
Thus the constrained donor search leaves an appreciable force on the frozen
metal. These are vacuum-MACE loads, not gradients of the solvent-composite score.
The force direction at an endpoint's own geometry is not a matched Ca−La
differential observable.

Testable hypothesis: allowing the ion and the selected intact donor groups to
accommodate jointly may remove an artificial fixed-metal constraint. Neither
lower native energy nor reduced force proves improved discrimination. Later
scoring, if authorized after preparation review, must expose both candidate
geometries to both metals with unchanged solvent accounting and retain the
earlier candidate pool.

## Frozen proposal definition

- The exact same four angular primitives selected at the original paired q0 by
  the completed angular protocol; do not select again at its final geometries.
- Add all three existing Cartesian metal-translation primitives as a common
  block, with no force-direction or label-selected axis.
- Both metals use the same seven IDs, with units explicitly `[Å, Å, Å, rad,
  rad, rad, rad]`. Translate the actual metal nucleus; preserve native endpoint
  charge, atom inventory, source preparation, cofactor, water and protonation.
- One original-q0 start per metal, eight endpoints. Minimize each metal's own
  native-vacuum MACE energy, never the Ca−La score. PQQ, waters and all other
  nonselected physical atoms remain fixed. No new starts or composite gradients.
- Existing analytic Kinematics and cap chain rules. Constrain every physical
  source heavy atom to final displacement ≤0.8 Å. This also constrains the
  metal's Cartesian displacement norm, not merely its three component bounds.
  Angular components additionally remain in ±0.8 rad; Cartesian components
  additionally remain in ±0.8 Å. The component bounds do not replace the sphere.
- Existing SLSQP, maxiter200, energy ftol1e−9 eV and unchanged source/overlap and
  final-state checks. Optimizer derivatives are eV/Å for translations and eV/rad
  for torsions; no vector is mislabeled as all-radian. The explicit mixed-unit
  optimizer parameterization is an algorithmic choice, not an energy correction.
- SLSQP may evaluate physically infeasible intermediate points within component
  bounds; retain their maximum extent and count. The 0.8 Å condition applies to
  final candidates. No clipping, added penalty energy or success substitution.

The earlier all-donor-plus-metal search with 0.20 Å/0.20 rad bounds already
worsened the tested native DFT gaps and reached its boundaries. This proposal
does not claim that adding metal motion is new. It isolates the three diagnosed
residual metal coordinates relative to the completed four-angle protocol, keeps
the four force-selected angles, and uses its already declared 0.8 Å final domain.
The previous negative result remains evidence against assuming benefit.

## Preparation checks and limits

Use the eight exact archived physical maps/forces. Verify projected gradients
against Cartesian work from geometry-only central differences (no displaced
MACE energies), including metal sign, units and cap mapping. Verify the full
seven-column constraint Jacobian by geometry-only differences on all real maps.
Check identical paired DOFs, source states, fixed nonselected atoms, the metal
sphere versus coordinate box, and explicit unavailable candidates before runs.

No new scientific energy or force calculation is part of preparation. No
curvature, entropy, populations, calibration, default change or guaranteed
stationarity is claimed. Root reviews this adapter and the completed 26-reference
angular continuation before deciding whether to run the eight starts. Any later
solvent calculation requires a separately explicit finite pool manifest.
