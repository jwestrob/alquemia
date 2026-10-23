# Correlated reference qualification: recommended next finite tier

21 October 2026. Design frozen before correlated calculations. No new native
calculation is released by this document. Root owns scheduling and resource
review. Keep the original PBE0 recipes/results unchanged.

## Question and priority

Does the opposite acetamide-versus-water metal contrast come from the interaction
model, or from an unreliable quantum comparator? Start with the smallest paired
water/amide systems, whose exact coordinates already have MACE and AMOEBA results.
Then test whether the answer survives explicit first-shell waters. Neither the
two PBE0 recipes nor the published La-only force field establishes La/Dy truth.

## A. Tractable closed-shell correlation pilot

Use the exact central water and acetamide dimers in prepared_v2/PANEL.json,
without optimization. Retain the archived La ECP46/Dy ECP55 and metal basis,
def2-TZVP light basis, effective singlets and separate physical La1/Dy6 metadata.

1. **Six canonical CCSD(T) water cells:** La and Dy, each complex, real ion with
   ghost water, and real water with ghost ion. This measures a counterpoise (CP)
   corrected interaction for each ion with only18 explicit electrons in the
   complete complex. It also establishes actual cost before larger cells.
2. **Twelve DLPNO-CCSD(T1)/TightPNO cells:** water and acetamide, La and Dy, the
   same three fragment roles. Water overlaps geometrically with the canonical
   pilot to measure the local-correlation approximation. Amide adds the paired
   donor-versus-water contrast. No D4 correction is added to coupled cluster.

Use an RHF reference, VeryTightSCF, no automatic restart/fallback, and explicit
NoFrozenCore for all electrons left outside the ECP. This preserves correlation
of the eight remaining metal electrons instead of relying on an implicit
element-specific frozen-core default. It does not recover the removed4f shell.
Record the actual correlated orbital count and the limitations of a valence
basis for remaining semicore/ligand-core correlation. Auxiliary bases for the
local method must be generated/pinned for both orbital bases before release.

The runner must prove that ghost centers retain basis functions but contribute
no nuclear charge, electrons or ECP potential. Real ion fragments have charge+3,
water/amide fragments0, all effective singlets. Verify native electron/basis
counts for every role; ghost metal handling is an explicit preflight item.
The intended interaction is E(complex, full basis)−E(ion, ghost ligand basis)
−E(ligand, ghost ion basis). Record the raw four-complex contrast as well as its
CP-corrected value; do not mix CP and uncorrected terms within a cycle.

## B. Basis and force qualification, after actual pilot review

For those same four dimers, a12-cell CP repeat with def2-QZVPP light atoms tests
light-basis sensitivity while retaining the metal basis. This is explicitly a
partial basis test, not a complete-basis limit. A separately sourced compatible
larger La/Dy metal basis is required before claiming full basis convergence;
do not invent exponents or change the ECP during that comparison. The acetate
extension must also test diffuse light-atom basis functions for its anionic
fragment and preserve separate mono/bi internal geometries.

Directional forces require new energies because DLPNO-CCSD(T1) has no analytic
gradient in this implementation. For water and amide, both metals, displace only
the metal along the original radial axis by±0.01 and±0.005Å. The paired two-step
CP design has48 displaced fragment-role energies, if all four centers remain
admitted. Compare the two finite-difference estimates, canonical water checks,
and original model/reference force projections. PNO thresholds can cause small
surface nonsmoothness; a failed step comparison stays failed rather than being
hidden by selecting whichever displacement agrees with a model.

## C. Admission, error allocation and state limits

Require native normal termination, correct physical/effective states and basis
inventory, finite energies, SCF and CC residual reporting, and a consistent
fragment/correlated-electron convention. Report T1/D1 or available equivalent
diagnostics, largest amplitudes, triples contribution and reference stability;
single-reference warnings prevent treating CCSD(T) as qualified truth.

For this exploratory tier, predeclare0.25kcal/mol as the desired stability of
the paired exchange contrast under local-correlation and available basis/CP
tests, and0.5kcal/mol/Å for radial-force step consistency. These are working
resolution targets, not measured error bars or replacements for experimental
uncertainty. Failure means the reference remains unresolved. Even passing them
does not close the missing metal-basis, frozen-f, spin-orbit or ensemble error.
Allocate the final reference error budget against the observable before any
held-out validation or production claim.

The explicit-f Dy PBE0 calculations test a different scalar representation, not
the complete f9 multiplet or spin-orbit problem. If its paired contrasts depend
materially on that representation, use an explicitly declared state-averaged
4f active-space/multireference+SOC study on the small paired ligand environments.
A candidate minimal4f active space is9electrons/7orbitals, but ligand-to-metal
charge-transfer character may require enlargement. Occupations, selected roots
and ligand-field/SOC shifts must be checked before release; a nominal sextet
determinant is not that qualification.

## Sources and interpretation

The author study supplies La structures and CP-corrected correlated interaction
curves; it supplies no paired Dy validation. [Wait et al. and retained SI](https://doi.org/10.1021/acs.jpcb.2c07237).
ORCA documents the iterative-triples keyword, single-reference limitations,
default frozen cores and lack of DLPNO analytic gradients. [ORCA6.1 coupled cluster](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/mdci.html).
Its ghost-center convention provides the required fixed-geometry fragment basis
tests. [ORCA6.1 counterpoise](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/counterpoise.html).

These are conditional electronic interaction checks. Hydrated fixtures, external
fields, sampled configurations, free-energy overlap and the four bound/water
quantum correction terms remain necessary for the complete selectivity goal.
