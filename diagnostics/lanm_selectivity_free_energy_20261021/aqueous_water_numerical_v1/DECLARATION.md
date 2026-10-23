# Frozen water candidate: transfer into the BIO parameter namespace

Declared 21 October 2026 before numerical execution. Root owns submission.

## Question

Does the exact frozen17605 paired water model retain its intended interactions
when translated into the existing BIO parameter file for the aqueous leg?
The known SI/BIO Urey–Bradley equilibrium-distance difference must remain visible.
This is a mapping check, not a new fit or a bulk-water physical qualification.

## Inputs and finite calculation

Use the independently admitted17605 water candidate, the exact17612 consumed
central aquo8 La and Dy geometries and native outputs, and the unchanged BIO
parameter file recorded in periodic_mapping_v1. The new aqueous_water_mapping_v1
manifest binds their full hashes and resolves metal types/classes900/901,
water O402/class73 and H403/class74. A common metal kinetic mass138.91 is a
recorded sampling convention; it does not enter these potential energies/forces.

Exactly two new nonperiodic25-atom cells: central aquo8 La and Dy. Each receives
one native CPU analyze and one analytic testgrad invocation, four programs total,
zero new finite differences. Reuse both existing SI results; do not recalculate
the independent orientation references or expose reserved donor outcomes.

Use the already pinned Tinker26.2 CPU executables, one numerical thread, mutual
polarization with polar-eps1e-8, exact frozen metal/water pair parameters and all
native BIO nonmetal terms. No periodic boundary, implicit solvent, optimization,
changed charge, fitting or dynamics. Standard_expand, one allocated CPU, default
memory/time. Preserve the native300s per-program limit of the existing adapter.

## Declared comparison

The only intended potential-energy difference is the water Urey–Bradley term:
force constant−7.60, equilibrium H–H distance1.5337Å in SI versus1.5537Å in BIO.
Read and pin actual global Urey coefficients and source formula before execution;
require unit1, cubic0 and quartic0 for this analytic reduction.

For each water H pair with distance r, calculate
deltaE = k[(r−1.5537)^2 − (r−1.5337)^2].
Its gradient difference on the first H is
2k(1.5337−1.5537)(r_H1−r_H2)/r, with the opposite on the second H.
The oxygen and metal contribution from this correction is zero. This correction
explains a known Hamiltonian difference; no BIO value is overwritten to force parity.

Require complete finite raw energies and all25 indexed gradient vectors, native
success with no failure/nonconvergence marker, exact input/source pins, and exact
chemical/type correspondence. After subtracting the calculated Urey difference,
require absolute total-energy difference≤1e-4kcal/mol and maximum absolute
Cartesian gradient difference≤1e-3kcal/mol/Å, matching the existing native-engine
parity limits. Require native Urey component difference to agree within1e-4 and
all other reported component differences≤1e-4. Both metals must pass; missing
or malformed output cannot pass. Retain uncorrected differences and correction.

## Interpretation and continuation

A pass establishes this finite namespace transfer, conditional on the frozen
small-system model. It does not establish flexible-water populations, periodic
electrostatics, salt association, hydration free energy, protein interactions or
La/Dy selectivity. Prepared full aqueous boxes remain unrun. Separate periodic
numerical and physical hydration tests are required before aqueous production;
the unfinished donor fit and independent donor/protein gates still govern the
bound leg. A failure is preserved and diagnosed; no automatic parameter change,
weakened gate, retry or promotion follows.
