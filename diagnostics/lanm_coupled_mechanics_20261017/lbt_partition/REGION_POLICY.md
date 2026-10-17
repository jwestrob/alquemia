# LBT local electronic region and native MM ledger

Module 3 of the activated parent plan. Preparation and real-fixture checks run
on Slurm. This adapter performs no electronic, solvent or MM energy evaluation.
It does not alter the existing whole-peptide model or its results.

## Region rule, frozen before new energies

Use the exact 256 physical atoms from preparation12469, including its conditional
free NH3+/COO− termini. No H normalization, missing-atom preparation or relaxation.
7CCO is La-conditioned; 7CCN is Lu-conditioned. All La/Dy coordinates are paired.

The baseline starts with complete residues3–14: every donor residue, the peptide
backbone between them, Asn6 and Asp13 are retained. Before constructing either
region, find direct hydrogen-bond candidates involving the ten declared donor
oxygens or Asn5ND2 in either original source. A donor is native N/O bearing H;
an acceptor is native O. The declared geometric screen is D–A≤3.5Å,
H–A≤2.6Å, D–H–A≥120°, excluding pairs within three native covalent edges.
This is a deterministic inclusion screen, not experimental H-bond assignment.
Include the complete residue of every external partner, then take the contiguous
sequence interval containing that union. Both sources use the same interval.

Each internal left boundary includes the preceding residue's C/O; each internal
right boundary includes the next residue's N/H. Cuts are therefore actual C–CA
and N–CA bonds, capped with H at1.09/1.01Å along the true omitted-anchor direction.
The two-anchor Jacobians are the existing `affordable_response.cap_jacobians`.
No legacy48 cap is reused, and no cap is an independent physical/MM particle.
The expanded region adds one complete residue at each available end and repeats
the same boundary rule. If baseline already covers the whole peptide, expansion
is explicitly unavailable; do not manufacture a distinct partition.

## Charge and interactions

The electronic charge is the sum of full-residue formal charges, independently
checked against native ff19SB residue totals, plus target+3. Partial formyl/amide
boundary fragments are neutral. Physical La multiplicity1 and Dy multiplicity6
are separate from GFN2 UHF0 / native frozen-f effective-singlet conventions.

Use the existing Mex residue-local shift: remove selected-QM and omitted-MM1 CA
charges from the embedding field; distribute the charge needed to preserve each
partial residue's original total equally over the two remaining heavy neighbors
of that CA. Record every removed/recipient ID and increment. Charges do not vary
with coordinates. The unshifted full ff19SB+formal3 solvent vector remains a
separate artifact. There is no cap population projection.

Reuse Mex's actual boundary ledger: omit wholly-QM bonded terms; retain wholly-MM
terms and every cross-boundary radial bond; omit cross angular/torsion terms
whose only MM support is represented by an actual cap, retaining other cross
terms. Record every native term and decision. Retain MM–MM LJ/Coulomb with native
exclusions/1–4 scaling, and real-QM–MM LJ with native exception replacements.
No classical QM–MM Coulomb, QM–QM FF nonbonded, cap FF term or C4 term is added.
Embedding uses all retained shifted point charges without native cross-region
electrostatic 1–4 scaling. Its boundary charge representation is the same Mex
approximation, not an exact reconstruction of the original full-FF Coulomb term.
The target uses the same explicit pure12–6 TIP3P La/Dy parameter file as Mex;
its transfer to the electronic boundary remains an approximation.

The eventual approximate energy is localMACE + GFN2embedded−GFN2vacuum + thisMM
ledger + one global fixedFF-ALPB term. EmbeddedDFT replaces both localMACE and
the GFN2 difference. Never add this ledger atop whole-peptide MACE.

## Preparation gates and next qualification

Both sources and both regions remain in the denominator. Verify exact source
order/coordinates, native topology/charges, all covalent cuts, donor/H-bond
coverage, formal/field/global charge closure, state parity, complete term and
exception accounting, and the source-derived cap/field chain rule. Two-step
finite differences of the map use1e−5/5e−6Å with maximum derivative error1e−6;
translation closure uses1e−12. These are geometry tests, not molecular forces.
Serialize the actual native MM systems and prove La/Dy share the same first
three component definitions. No Context is created in preparation.

The ledger is not energy/force qualified by these checks. Root must bind the
prepared output, then qualify native MM and complete electronic+MM+solvent
derivatives and same-geometry region sensitivity before optimization. This
module authorizes no new molecular calls, classification, parameter fit or
production change.
