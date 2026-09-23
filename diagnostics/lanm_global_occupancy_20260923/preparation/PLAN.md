# Full-chain La/Dy occupancy preparation

Jacob authorized a small within-lanthanide pilot with two- and four-ion states,
starting with Hans. Parent assigned exactly Hans 8DQ2 A, Hans 8FNR A and Mex
8FNS A; for each prepare EF1+EF2, EF2+EF3 and EF1–4, with La and Dy endpoints.
This script performs preparation/audits only: zero energy, optimization, new
protonation, folds or Slurm submissions. Existing pH5 standard-residue H inputs
are pinned and reused. No reserved SpyCI assay outcomes are consulted.

Retain the complete observed chain A (110 Hans residues or 105 Mex residues),
all its author-chain crystal waters, and no other chains/heteromolecules.
No missing heavy atoms, loops, terminal atoms or peptide caps are synthesized.
The existing ff19SB topology audit supplies covalent connectivity and integral
protein charge. Its established radial protein-H projection fixes bond lengths
without changing H inventory, attachment or direction; water H directions are
likewise retained with TIP3P O-H length 0.9572 Å. All heavy coordinates must match
the archived selected experimental source exactly. Record every H adjustment.
No angle/coordinate relaxation is used. Source water inventory remains identical
across occupancies/metals within each source; it need not match between crystals.

EF1/2/3/4 map to actual source residues A201/202/203/204. 8DQ2 has La/La/La/Na;
8FNR has four deposited Dy; 8FNS has four deposited Nd. Remove the nonselected
ions for two-ion hypotheses; substitute every selected ion with the declared
endpoint element. Thus four-Ln 8DQ2 replaces source Na at EF4. These conditional
monomer states do not assert natural four-Ln occupancy or Hans dimer equilibrium.

Total charge = audited protein formal charge +3*n, with neutral retained waters.
Physical La multiplicity1 and maximum-spin Dy multiplicity1+5*n are declared;
other inter-ion couplings are untested. Native GFN2 uses an effective singlet,
with exported three-electron 5d/6s/6p lanthanide parameters and electron parity
checked independently. No physical f electrons are passed to the f-in-core model.

Deliver nine states/eighteen paired XYZ files, complete atom mapping and covalent
bond graph, source/protonation/parameter pins, completeness and geometry audits,
and explicit failures. Protein/water atoms precede occupied ions; metal indices
are supplied, never inferred from ordering. No energetic labels or bands apply.
