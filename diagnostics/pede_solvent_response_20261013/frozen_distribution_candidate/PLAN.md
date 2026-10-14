# Frozen-distribution candidate: predictive test, not a fitted response factor

Declared14October2026UTC after the real eight-cell DFT check and native ALPB
fixed-charge replay. Jacob's contained solvent-correction development remains
authorized. Production, thresholds, old failures and all other sessions stay
unchanged. This document fixes the model and evaluation sequence. Execution
waits for the exact finite source/cell manifest and existing calibration-rule
pins now being prepared; no candidate scores have been inspected.

## Question and model

Does omitting the cheap model's additional solute electronic redistribution
improve reference decisions/source robustness without losing reference fidelity?
The replay located18.7467 of21.2566kcal/mol source sensitivity in redistribution
plus intrinsic compensation. That component is not proven erroneous: independent
DFT retains10.0342kcal/mol and distributed charge response. This candidate is a
controlled approximation to test, not a claim that all redistribution is false.

For every metal and geometry in each existing compatible shared pool:

`E_candidate(M,q) = E_OMOL,vac(M,q) + G_nativeALPB(q_GFN2,vac(M,q); geometry q)`.

G includes native electrostatics, charge-dependent hydrogen bonding, SASA and
reference terms. The native reaction field responds to fixed solute charges;
the solute density does not re-equilibrate. This is a perturbative mixed-model
descriptor, not exact replacement of self-consistent solvation or binding free
energy. No interpolation, fitted response factor or automatic inherited offset.
It is NOT the old self-consistent U with a forgotten term: it is an explicitly
different model with its own identity and calibration.

Reuse only exact compatible archived MACE energies, native GFN2 vacuum/water
charges, coordinates and states. No new SCC, DFT, gradients, optimization or
folding in this phase. Missing artifacts stay missing. No charge reassignment,
global neutralization, altered proton/water inventories or parameter tuning.

## Shared candidate selection and calibration

Each source keeps its complete origin/Ca-proposal/La-proposal pool and any other
candidates already admitted under the established source protocol. Both metals
must be scored on every required compatible geometry, with their own charges.
Do not retain whichever minima the old model selected. Keep mathematical row
minima and the existing numerical selection/tie policy explicit. A failed
required cell makes that pooled score unavailable; never substitute a successful
baseline field into the candidate. Do not mix absolute energies across folds
with different fragments or physical states.

First use the designated25canonical training inputs and the exact existing
calibration RULE, pinned by the inventory. The rule can yield a new scale, but
must not change in response to outcomes. No successful calibration is claimed
if required canonical inputs are missing or its separation/admission gate fails.
Stop dependent transfer scoring if that gate fails; report all denominator cells.

Only after the new calibration passes, freeze its result before structural
transfer. The finite transfer population must include the existing corrected
reference cohorts, full125Ca-source denominator and all15two-Ca folds, with exact
deduplication and recorded recovery identities. Inventory will pin the available
La-source cohort rather than infer it from old prose. These repeatedly examined
structures are development evidence. No new reserved outcomes or PLM labels.

Compare against the strongest compatible corrected baseline using unchanged
aggregation/grouping rules, reporting correct/wrong/inconclusive/unavailable,
common coverage, individual source contrasts, selected geometries and group
spread. Report old-band transfer separately where meaningful; do not refit on
transfer outcomes, choose favorable source subsets or call abstention accuracy.
No whole-library/PLM production rescore or promotion follows automatically.

## Numerical and execution gates

Reuse the successful v3 xTB6.7.1 build and native model, with an explicitly
mapped variable atom count. Compile only the thin wrapper against pinned native
libraries/modules onSlurm. No new solver implementation or environment install.
For each endpoint, replay its saved water charges against all native components
within the unchanged precision-derived bounds before interpreting its vacuum-
charge evaluation. Preserve cell-level failures and physical state/charge checks.
Required-cell completeness applies to every pool and calibration.

Use actual saved charge precision; never invent shell populations or silently
relax bounds. Preserve current native temperature, radius/reference conventions
and corrected original parameter lineage. Fortran input conversion must use the
actual mctc_io XYZ-reader convention qualified in v3. Every cache identity includes
state, geometry, charge source, native implementation and candidate protocol.

All build/tests/evaluations useSlurm under current shared allocations; root owns
submission and durable completion wakes. Measure actual preparation/build/replay
and collection costs separately from reused expensive charge/MACE generation.
The planned run is a finite archived-artifact evaluation, not permission for
new quantum endpoints or indefinite parameter searches.

Curie owns exact population/cell/calibration-rule inventory. Heisenberg owns the
thin native replay adapter. Root freezes the integrated manifest, runs the stages,
interprets results and updates shared status. The first stage may legitimately
reject this candidate; do not tune the model or admission rule to force separation
of the inspected labels.
