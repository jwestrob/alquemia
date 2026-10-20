# One directional consistency check of the whole-peptide PBE0 reference

Declared 20 October 2026 after native15545 and saved collector15838. The
model-converged step725 structure retains PBE0-plus-fixed-ALPB force RMS3.587757
kcal/mol/Å versus model0.029611. Reference and model support the step92-to725
energy lowering but differ by1.228612kcal/mol in work. These results establish
disagreement; they do not qualify the reference derivative or identify which
electronic Hamiltonian is more accurate.

## Question and frozen selection

Does the admitted whole256 PBE0 vacuum force at the exact step725 structure
agree with a centered energy derivative along the direction of the full
model/reference force discrepancy?

Use native RESULT `a53125b98e597ea0473e06d58316ac2b32c914e9647da595f163cf0aa76a8c19`,
collection RESULT `407bdf88e42c43cb2ae48ed1493cbbc008f5c0eb957628787a6aa1c1b795bb0d`,
and the new_La force report
`e6c0d039500c9deb87a53f31f974dfdc741c1cd801e481b51dc3316d100b7a7c`.
Reuse the exact256-atom coordinates, IDs, whole charge−2, La physical/effective
singlet and source chemistry. Retain all atoms; no alignment, atom selection,
fitting, rigid-body projection or geometry optimization.

Define `u = (F_MACE,vac − F_PBE0,vac) / ||F_MACE,vac − F_PBE0,vac||` using the
Euclidean norm over all256×3 components. The shared solvent cancels from this
discrepancy. Generate exactly two geometries `x±h*u`, with **h=0.01Å in the full
Cartesian vector norm**. No single atom can move more than h. Slurm preparation
must verify finite vectors, exact saved identities and units, norm/sign, source
state and unchanged atom inventory before freezing the two input decks.

## Native evaluations

Exactly two independent whole256 vacuum singlepoint energies, plus/minus. Reuse
the15545 native input and its pinned executable/La basis/ECP46/auxiliary basis,
PBE0-D4, def2-TZVP, def2/J, RIJCOSX, DefGrid3, VeryTightSCF, NoAutostart and
PModel guess. The only deck changes are removal of EnGrad to request a
singlepoint and the explicitly declared displaced XYZ coordinates. Retain
actual native convergence semantics and all printed residual/gap warnings;
do not replace native admission with a favorable energy-only inference.

Use12MPI per endpoint, standard_expand, scheduler default memory/time, and
allocation-derived MaxCore. Preserve the existing runtime renderer, receipts
and hashes. No central rerun, wavefunction seeding, SCF/grid adjustment, new
MACE/ALPB/GFN/MM call or other-metal calculation. Arm actual completion and
alert-only health wakes. Each endpoint runs once after source review and the
oneCPU zero-call preparation has passed. Retain unavailable endpoints.

## Comparison and decision

The saved-only oneCPU collector subtracts the final printed native energy
strings using decimal arithmetic before unit conversion. Compare
`F_FD = −[E(x+h*u)−E(x−h*u)]/(2h)` with the saved central `F_PBE0,vac dot u`.
Record both values, their signed difference, both source energies, all
admission/SCF diagnostics and actual costs. Reproduce the central gradient's
reported units/order and the full model/reference directional projections.

Freeze the numerical consistency rule before executing either displaced SCF:
`abs(F_FD−F_analytic) <= 0.005 + 0.0001*max(abs(F_FD),abs(F_analytic))`, in
kcal/mol/Å. This reuses the earlier project's directional reporting rule; it
does not relax a previous gate or define physical accuracy. A pass establishes
only consistency for this sampled direction and step. A failure or unavailable
endpoint stops this finite test without automatic retries, half-step sweeps,
new grids or tolerance changes. Neither outcome validates all force components,
electronic state, an equilibrium minimum, solvent, metal preference or transfer.

The prior Mex external-charge derivative failure is motivation to verify a
reference, not evidence that this whole-system nuclear derivative is defective.
MACE-OMOL's training Hamiltonian differs from this PBE0-D4 comparator; preserve
that distinction before attributing the residual discrepancy to learned forces.
The test is separate from the closed Mex branch and does not reopen it.

## Numerical context

ORCA's documentation distinguishes the COSX grids used during SCF and final
energy/gradient evaluation and notes grid noise in numerical derivatives.
Here those settings stay fixed so the test assesses the existing reference.
See [ORCA6.1.1 numerical integration](https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/numericalintegration.html)
and [numerical gradients](https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/numericalgradients.html).
One sampled derivative cannot by itself separate SCF, grid, finite-step and
electronic-state causes if the check fails.
