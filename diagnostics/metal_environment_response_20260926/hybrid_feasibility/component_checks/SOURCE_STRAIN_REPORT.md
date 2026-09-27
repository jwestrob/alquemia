# Large classical components: archived source strain, not a found unit defect

**The archived 1H4I preparation contains severe reconstructed-tail clashes and
systematically overlong hydrogen bonds.** Read-only checks found no mismatch in
retained native force constants, indices, LJ parameters or exceptions. These
are concrete source defects for mechanical response; this is more specific than
calling a large total energy ordinary whole-protein strain. It does not uniquely
decompose every kcal/mol in the combined bonded force group.

No new potential evaluation, optimization or molecular job was performed.
Root completed the frozen derivative matrix unchanged:32/32 configurations and
all200 classical checks pass (`FINAL_COLLECTION_1219497.json`). Correct
derivatives do not make these coordinates a mechanically equilibrated source.

## Actual scout and exact parameter audit

Ca_A job1219495 returned bonded55,319.427213, MM LJ11,604.178599,
MM Coulomb−14,673.884719 and cross LJ−57.146019 kcal/mol. These are existing
scout outputs, not fresh calculations. Its maximum saved MM-LJ force is
94,153.08 kcal/mol/Å on the rebuilt A596 O atom.

Comparing serialized native parameter objects directly verifies:

- All9,225 retained bonds,16,655angles,23,666torsions and597CMAP supports and
  coefficients are identical to the selected original terms. All16CMAP tables
  are exactly preserved.
- All9,090MM LJ particle parameters and49,765MM exception LJ parameters are
  exactly native, including exception particle indices.
- Native lengths are nm and input coordinates convert Å→nm once. Source bond
  geometry is on the expected Å scale; no bond exceeds its native equilibrium
  length by more than0.239Å. This argues against a gross tenfold distance error.
  The code uses quantity conversion to kcal/mol and kcal/mol/Å once.

This rules out the tested mapping/parameter/unit failure mechanisms, not every
possible implementation defect. The completed finite-difference and rigid checks independently support the
classical implementation, despite the strained input.

## 1. Nonexperimental tail is severely clashing

Original `workspaces/mxaf_qm/1H4I.pdb` REMARK465 lists chainA residues596–599
(SER,ALA,ALA,LYS) as experimentally missing. The old
`1H4I_protonated.pdb` contains these residues with26 added heavy atoms including
the terminal OXT. Its header dates preparation to7May2026/OpenMM8.1.2. We have
not reconstructed the exact old repair command from that header.

All4,620 common chainA heavy atoms, including nonprotein entries, have exactly
identical deposited and protonated coordinates. The problematic tail atoms are
added coordinates, not experimentally resolved atoms displaced by this ledger.

| Actual source pair | Separation | Interpretation |
|---|---:|---|
| A595 CE – rebuilt A596 O |1.406364Å| Nonbonded C/O clash; no native exception |
| A595 CG – rebuilt A596 CA |1.843958Å| Nonbonded C/C clash; no native exception |
| rebuilt A599 CB – A599 O |1.665282Å| Compressed native1–4 pair; its scaled LJ exception is retained |

The exact same separations are present in the archived protonated PDB and
prepared coordinates. The first pair's atoms carry the two largest saved MM-LJ
forces (94,153 and92,202 kcal/mol/Å). They lie about56Å from the metal; other
tail clashes lie54–65Å away. This is evidence of a remote preparation artifact,
not metal-pocket strain to be interpreted as selectivity. No pair energy was
newly evaluated, so an exact percentage of total LJ energy is not asserted.

## 2. Hydrogen geometry is systematically incompatible with ff19SB

Retained source bonds show:

| Bond type | Count | Median actual length | Median excess over its native equilibrium |
|---|---:|---:|---:|
| C–H |3,451|1.192150Å|+0.103088Å|
| N–H |914|1.185522Å|+0.175522Å|
| O–H |92|1.185788Å|+0.225788Å|

Every observed C–H,N–H andO–H bond is longer than its native equilibrium.
For comparison, median heavy C–C/C–N/C–O residuals are about−0.00003/−0.00132/
+0.00248Å. This is an element-dependent hydrogen-preparation mismatch rather
than a global coordinate-unit error. Angles involving H also have median
absolute equilibrium departure6.46° and maximum65.43°, versus1.04° median for
heavy-only angles. The largest saved bonded forces are at amide/amine nitrogens
across the protein, consistent with widespread hydrogen-coordinate strain.

We did not calculate bond/angle energy sums here, so the entire55,319kcal/mol
cannot be assigned exclusively to these hydrogen defects. Nonetheless, using
this source as a mechanically equilibrated protein would be indefensible.
The current `protonate_cif.py` does not reconstruct missing residues by default;
that current policy must not be retroactively credited to this May artifact.

## Fixed local contrasts versus relaxation

The completed matrix finds **exactly zero classical Ca/La differential A→B
response** in every component for this fixed Thr159 hydroxyl perturbation.
Ca and La share the same MM model and A/B movement; metal-position and other
metal-LJ geometry are unchanged. Thus the large common classical totals do not
create or obscure a classical metal-dependent response in this test. They are
not evidence that the existing electronic double difference must be wrong.

That cancellation does not mean the source defects are irrelevant to the
electronic calculation: the distorted hydrogen inventory is the actual embedding
geometry, including the locally perturbed hydroxyl. The unchanged remote tail
can also contribute an electrostatic field. Quantifying either effect needs a
separately declared repaired-source comparison, not an assertion from distance
or energy magnitude alone. No such molecular comparison ran here.

For relaxation the problem is immediate: ~56Å-distant terminal clashes and
thousands of strained H coordinates would produce major motions unrelated to
local metal accommodation. Allowing this preparation to relax would confound
source repair with metal-dependent response. A corrected force implementation
can faithfully drive an inappropriate starting structure.

## Implication for the next solvent-consistent test

1. Preserve this scout and its passed classical derivative checks as an immutable
   test of the declared input; do not reinterpret that as full hybrid qualification.
2. Prepare a new common physical source before investigating scaffold response:
   retain experimentally supported heavy geometry, repair hydrogen coordinates
   under the declared proton inventory, and resolve the missing tail explicitly.
   Options include a documented resolved-chain construct with explicit terminal
   chemistry, or independently justified missing-segment preparation. Do not
   silently delete the tail, change total charge, or treat generated coordinates
   as observations. Whichever construction is adopted must be shared by metals.
3. Introduce solvent in the **same force-generating Hamiltonian** before any
   protein relaxation. Solvent addition alone will not undo these hydrogen and
   missing-tail defects, and an old GB score cannot supply consistent forces.
4. Validate the repaired source and additive coupling, then assess metal-dependent
   physical response. A large source-strain relaxation is not itself a selectivity
   improvement. Hydration-derived metal LJ remains unqualified for QM/MM.

No optimization of the present dry model is recommended or implemented.

## Reproducible read-only evidence

- `SOURCE_DIAGNOSIS.json`: pinned source/scout, actual geometry statistics,
  original missing-residue remarks, pair separations and saved force extrema.
- `PARAMETER_MAPPING_AUDIT.json`: exact native-to-scout serialized parameter checks.
- `DIAGNOSE_SOURCE.py` and `DIAGNOSE_PARAMETERS.py`: reproduce those observations
  without OpenMM Contexts or potential evaluations. Their existing output files
  are immutable; use a new destination/version for any later regeneration.
