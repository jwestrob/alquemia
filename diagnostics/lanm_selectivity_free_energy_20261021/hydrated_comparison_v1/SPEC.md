# Saved comparison of hydrated donor competition

Question: how does including first-shell water change the paired La/Dy donor
response across MACE, transferred AMOEBA, and the frozen-f PBE0 reference?
Use the declared 18 mixed fixtures, six central native cells, and two previously
computed aquo8 controls. No new molecular calls, fitting, or geometry changes.

Verify the original and mixed native panel provenance separately, every input
geometry and chemical state, result hashes and native admission checks. Report
missing or failed references explicitly; never fill them with another method.
All arithmetic uses one local CPU within root's seven-thread share.

For each donor and radial position report E(Dy,mixed)-E(Dy,aquo8) minus
E(La,mixed)-E(La,aquo8). Free ligands and displaced waters cancel in this balanced
contrast. Negative means conditional Dy stabilization relative to the specified
aquo8 configuration; it is not a binding free energy or equilibrium preference.

Compare full forces and the derivative along rigid ligand translation. For the
latter use minus the sum of forces on ALL ligand atoms dotted with the declared
translation axis. The metal-only force projection is not this derivative.
Report metal-difference force errors separately from individual-metal errors.

The isolated and hydrated fixtures also differ in surrounding geometry. Their
comparison tests model behavior across these contexts, not a unique additive
hydration contribution. Native SCF admission is not strict residual, electronic
state, basis, spin-orbit or reference-accuracy qualification. No new acceptance
threshold is selected from these outcomes.
