# Frozen final water counterpoise comparison

Declared before inspecting any new DFT counterpoise output. Saved-only collection;
no molecular calls, new solver settings or modified native admission.

## Admission and bookkeeping

Require all five DFT fragment results, both original PBE0-D4 complexes, all six
canonical CCSD(T) results and the complete canonical collection. Verify exact
source/panel/input/output pins, per-cell inventory and executed input decks.
Replay the frozen native parsers on actual outputs; every originally admitted
field must match. Recompute the canonical collection and use the existing
DFT collection arithmetic. Retain all terms in Hartree and the final units.

PBE0 electronic energies are separated from real-atom D4 dispersion. No D4
evaluation on ghost atoms is admitted. The real water correction enters both
individual interactions, and isolated-ion dispersion is the recorded exact
empty-pair/triple sum. Retain electronic PBE0 CP and PBE0-D4 CP separately;
canonical CCSD(T) receives no dispersion addition. Report Dy-minus-La contrasts
and both individual-metal differences for the exact same fixed coordinates,
orbital bases and frozen-f effective-core representation.

## Finite decision

Primary comparison is the paired PBE0-D4 CP contrast minus canonical CCSD(T) CP.
Agreement means its absolute difference is at most0.25kcal/mol, adopting the
existing correlated-qualification exploratory resolution for this diagnostic.
Otherwise PBE0-D4 is not interchangeable with canonical CC at that resolution
on this fixture. Report the corresponding electronic-only PBE0 difference
descriptively; do not choose a Hamiltonian after observing which agrees better.
Individual-metal discrepancies remain visible even if their paired difference
cancels. This is method agreement at one basis/geometry, not physical accuracy,
an experimental error bar, donor selectivity, forces, or a binding free energy.

## SCF interpretation

Retain raw numeric residuals, native stopping signal and declared convergence
mode for every cell. Do not impose an all-six-printed-threshold requirement.
The saved mode review shows Total+1el-Energy/SOSCF and stale post-SOSCF DIIS
fields in the prior14 endpoints. New outputs are reported using their actual
signals; that observation is not blindly assigned to future outputs. Numerical
precision is a distinct future test, with the bounded proposal preserved in
the pinned SCF interpretation report. Basis, frozen-f, spin-orbit and ensemble
uncertainties remain open regardless of this comparison's decision.
