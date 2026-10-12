# Component-qualified factorial v2 — frozen before the remaining matrix

Root authorized this separate development protocol on 12 October2026 after the
four-input scout and four-call s9 diagnostic. It preserves v1's failed absolute
energy gate. No geometry, Hamiltonian, population state or classification band
changes. The actual composite is OMOL(vac) + GFN2(ALPB) − GFN2(vac), so its GFN2
input is a **same-metal, same-coordinate solvent-transfer component**.

## What the diagnostic established

Setting only the non-self-consistent three-body dispersion factor s9 to zero in
both independent implementations reduces their Ca/La vacuum total-energy
residuals to 1.55e-7/1.48e-7Eh, below the existing1e-6 gate. Their original
residuals were6.47e-5/6.73e-5Eh; electron populations agree within3.1e-8e.
This identifies that component to the gate's precision. The exact upstream
implementations use different three-body cutoffs (xTB40Bohr; tblite25Bohr).
No cutoff patch was executed, so cutoff choice alone is not independently proved
to explain every last digit. The term is geometry-dependent, not an atomic
energy zero or a metal-specific constant.

In xTB6.7.1 `d4_gradient`, the ATM term uses `zerovec` assembled from zero-charge
reference weighting; the corresponding weights do not depend on the solved
atomic charges. In tblite0.7.0 `get_dispersion_nonsc`, atomic charges are explicitly
zero before constructing C6 and evaluating the three-body term. Both use
coordinates, elements and fixed dispersion parameters, without an ALPB input.
Consequently, within each backend and a pinned geometry/state, write
E(medium)=E_rest(medium)+E_ATM(q): E(ALPB)−E(vac) cancels E_ATM(q) exactly.
We compute this difference from **actual full endpoints**, not fitted or inferred
dispersion offsets. We do not compare or minimize GFN2 totals across geometries.

## Frozen qualification and scope

For each metal/context/geometry, require exact identical XYZ hashes, charge,
unpaired count, temperature and parameters between vacuum and water. Both
backends must converge with the existing accuracy0.0001/max500 policy, charge
closure≤1e-5e, valid finite per-atom charges and declared electron occupation.
Every input's maximum independent all-atom charge discrepancy must be≤1e-5e.
For each metal, require
abs((E_ALPB−E_vac)_xTB−(E_ALPB−E_vac)_tblite)≤1e-6Eh.
All available exact archived primary-energy checks must also pass≤1e-6Eh.
Report every absolute-energy residual, including its original failed v1 gate.
No retry or gate relaxation; a failed pair's component stays unavailable.

Reuse the four actual v1 scout inputs and run exactly the remaining76 real inputs
on each backend,152 new calls maximum. The s9 diagnostic is additional and kept
separate (four calls). Two one-thread workers maximum; no MACE/ORCA/optimization.
Unrounded primary standalone xTB energies alone feed the composite. tblite is
the independent component check. Candidate membership, selection tolerance and
class reference are owned by root and are unchanged by this solver protocol.

Manifest and receipts: `workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2/`.
Sources: https://github.com/grimme-lab/xtb/blob/v6.7.1/src/disp/dftd4.F90 and
https://github.com/tblite/tblite/blob/v0.7.0/src/tblite/disp/d4.f90 .
