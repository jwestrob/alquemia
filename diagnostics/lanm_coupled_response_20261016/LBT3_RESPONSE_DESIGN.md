# Deterministic force-selected whole-peptide response

Declared16October2026 before inspecting the new composite force rankings.
This extends LBT3_SOLVENT_PLAN.md; no thresholds, biological labels, source
selection or original chemical state are fitted by this design.

Read both completed real-source64-mode catalogues from12509 and the admitted
four whole-origin composite energies/forces from12510 when available. If a
required force component is unavailable, the dependent profile is unavailable;
no vacuum-only replacement is called a composite result. NativeMACE12500
provides proposal-Hamiltonian force diagnostics separately.

For each source, flatten the analytic physical Jacobian column J_i over all256
real atoms and normalize it by its Euclidean displacement norm. Project each
metal's force independently and its difference onto this normalized direction.
Report nativeMACE, solvent and composite projections in consistent units. A
torsion angle and translation distance are never directly compared without
physical-displacement normalization. Large common loads remain visible.

Select the same small set of named coordinates for both sources and both metals:
retain all three target translations, then at most four torsions ranked by the
maximum of |load_La|,|load_Dy|,|load_Dy−load_La| across BOTH sources. Break exact
ties by the saved mode order. Greedy admission requires a residual normalized
physical-displacement norm>1e−5 against admitted columns in BOTH source metrics.
This prevents redundant directions without inventing stiffness or eigenvalue
corrections. Save the complete rankings and rejected directions. The basis is
selected by forces; any later search must minimize each metal's own energy,
not optimize their difference.

## Small derivative check before optimization

Use7CCO as the fixed first source. Nominate its nonmetal torsion with largest
absolute normalized differential COMPOSITE force, with ties by mode order.
This selects a real response direction without using an affinity label or
favorable energy. Whether it enters the small common search basis is separate.

Generate deterministic ±0.001 and±0.0005angstrom maximum-atom-displacement
rotations along that one actual graph mode, for both metals at identical
coordinates. Keep every bond, protected group and proton/solvent inventory;
check the existing kinematic Jacobian and no newly severe contacts relative to
the source. These are explicitly modeled numerical tests, not observations or
thermal samples. Keep origin coordinates and old failures immutable.

The finite derivative matrix is eight MACE energy calls and16standaloneGFN2
energy calls (eight geometries/states × two media), no new gradients or DFT.
It reuses the four admitted origin analytic forces. Preserve all backend
components and compare their derivatives individually and in the complete
energy. Fixed composite residual tolerance is0.2kcal/mol/angstrom plus5%of the
derivative magnitude, as declared in the parent plan. Report both step results,
their difference and residuals. Do not enlarge tolerances or reinterpret a
missing component as a pass. Derivative qualification is not biological accuracy.

This file authorizes staging the geometry/profile adapter now; actual molecular
submission requires the complete parent solvent result and physical preparation
checks. A separate bounded accommodation design follows the observed forces and
numerical check. No unbounded optimizer, global protein campaign, dimer campaign,
whole-LanM vacuum retry or PQQ change follows from this preparation.
