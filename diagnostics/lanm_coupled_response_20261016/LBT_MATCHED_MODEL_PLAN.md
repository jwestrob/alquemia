# Fixed-FF solvent matched Mex–LBT3 static comparison

Declared before any new call, 16 October 2026. Root owns staging and submission.
This is a finite consumed-source diagnostic; it opens no new library or reserved
outcomes. Source writing and review do not constitute execution.

## Question and denominator

Does the conditional Mex–LBT3 balanced La/Dy exchange retain its direction when
both systems use the same additive energy rule and fixed native FF-charge ALPB
functional? Retain **all four** LBT3 cells: model1 7CCO and 7CCN, each La and Dy.
7CCO was deposited with La and 7CCN with Lu; both receive the already prepared
La/Dy substitutions. No preferred source is selected after scoring.

Reuse the exact four whole-256 native vacuum MACE records from Slurm12500 and
the two admitted Mex static assemblies from Slurm12685. Make exactly four new
256-atom native fixed-charge ALPB evaluations, one per LBT cell, sequentially on
one CPU. There are no new MACE, GFN SCC, MM, DFT, gradients used for response,
geometry preparation, optimization, or explicit-water calls. The native wrapper
also writes its analytic gradient; archive it without extending qualification.
No retries or replacement cells. A failure retains its declared slot and withholds
dependent contrasts. Root records actual scheduler cost separately.

## Hamiltonian and partition

Use E = MACE_vacuum(Q) + [GFN2_embedded(Q; field) − GFN2_isolated(Q)]
+ retained_MM_ledger + native_ALPB(full physical geometry, fixed native FF q).

For LBT3, Q is the entire peptide plus target (256 atoms). The external field,
MM region and caps are empty, so the field-difference and retained-MM terms are
partition-defined zeros. They are not evaluated zeros. Do not add a whole-peptide
FF vacuum energy to whole-peptide MACE. Mex retains its exact admitted capped209
local Q, shifted field, MM ledger, full2049 cavity, three Nd spectators and source
waters. The construction is matched; different Q partitions are not asserted
to be numerically invariant or equally accurate.

LBT charge is −2 and physical MACE multiplicity is La1/Dy6. Its fixed solvent
vector joins the existing255 PEPTIDE_ATOMS.json native ff19SB charge records to
PHYSICAL.json by physical ID, then appends formal +3 for the actual target.
These charges are unmodified, not fitted, normalized or derived from GFN/MACE.
The existing preparation hypothesis is free NH3+/COO− with Asp3/7/13 and
Glu11/14 deprotonated: peptide−5, holo−2. The primary SI terminal/−4 ambiguity
remains unresolved; this comparison is conditional on the declared preparation.

Use the exact Mex native Fortran and build inputs reached through12566, changing
only its atom-count admission from2049 to256. Water ALPB, internal GFN2-xTB/ALPB
parameters, p16 kernel, gsolv state,298.15K,230angular points and zero ionic
strength/radius remain unchanged. No SCC calculator is instantiated. Preserve
getEnergyParts and gradient initialization. dq=0 is exact-declared-FF-charge
bookkeeping, not zero physical/model uncertainty. Serialize q and coordinates
with17significant digits and require float round trips. Old self-consistent LBT
solvent offsets and old frozen-GFN-charge values are excluded.

## Admission before the first call

Pin implementation, this plan, preparation12469, MACE12500, Mex12685 and their
exact configs, maps, raw component artifacts and qualified native build lineage.
Verify both implementation paths if a workspace snapshot aliases an original
source by SHA. Require actual parent success and exact declared denominators.

For each LBT source check256 unique physical IDs,255 native FF records, target
index255, sequence/terminal additions, no caps/waters/spectators, preserved source
H/heavy atoms and source mappings. Check native XML particle count and every
NonbondedForce q against the saved FF export; compare charge/element/residue/order
and coordinates against PHYSICAL and exact MACE XYZ/state. The existing recorded
peptide bonds and source pins remain unchanged. Charge closure tolerances are
1e−6e for255 FF charges and1e−5e for native CHARGE_SUM; coordinate replay tolerance
is1e−12Å for serialization roundoff. These are input checks, not fitted gates.

Read saved Mex MACE/GFN/MM/ALPB records and verify their raw artifact/map pins,
source/state binding and assembled energy identity. Reconstruct saved energies
with the original operation order, requiring exact assembly equality; preserve
the admitted MM internal component-sum check at1e−8kcal/mol. Do not rerun its
scientific parser or replace its accepted mapped forces. Confirm shared MACE
checkpoint/head/float64 settings and the identical native solvent functional.
All four LBT input admissions and saved Mex replays precede compilation/calls.

## Exchange, components and common source pool

For each source s report

R_s = (E_Mex,Dy − E_Mex,La) − (E_LBT,s,Dy − E_LBT,s,La).

This is the energy of Mex–La(+2)+LBT–Dy(−2) → Mex–Dy(+2)+LBT–La(−2).
One La, one Dy, total charge0, Mex spectators/waters and each source membership
are conserved. Positive R favors the left side, conditional relative La in Mex
and Dy in LBT. Calculate component differences in native units before conversion
(23.06054783061903kcal/mol per eV;627.509474kcal/mol per Hartree), then sum with
math.fsum. Report all four additive components and the four ALPB subcomponents.

A common two-source LBT pool is reported only if ordered physical identity,
sequence, terminal state, exact FF q vector, native XML topology, covalent graph,
whole-Q partition and electronic-state/model recipe match across both sources.
Coordinates may differ. The exact same two-source menu is then available to
each metal; require all four cells. Each metal selects its own minimum composite
energy, fixed7CCO-before7CCN order for exact ties, with all ties disclosed. Report
the resulting balanced exchange and selected-state component decomposition.
Retain both source-specific R values regardless of pool selection. This pool is
a finite common menu, not an equilibrium ensemble or population estimate.

No protein Kd, experimentally matched affinity, successful classifier, independent
blind-test or bulk-solvent accuracy claim follows. Mex is a fixed EF3 substitution
conditional on three Nd spectators; LBT3 is a monomeric peptide tag. Existing
source, hydration, f-electron/spin and partition limits remain.
