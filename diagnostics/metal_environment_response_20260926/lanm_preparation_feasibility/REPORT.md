# LanM EF3: reuse normalized sources, prepare a new local embedding

**The needed hydrogen repair already exists for all three consumed LanM sources.**
Reuse the whole-chain `prepared_v1/sources/*/atoms.json` coordinates and graphs,
not the original long-H protonated PDB coordinates or the failed relaxed proposals.
No selected-chain rebuilt-tail defect was found. Complete EF3 donor chemistry and
specific second-shell support can be represented in proposed195/198/209-atom
regions. These are source-membership proposals, not prepared or scored QM inputs.

No new molecular calls, protonation, coordinate repair, optimization, or reserved
outcome access occurred. Historical failed whole-chain jobs remain untouched.

## Actual source quality

| Source | Selected chain | Protein atoms | Full chain waters | Largest old-core H change already applied |
|---|---|---:|---:|---:|
| Hans8DQ2 | A24–133,110residues |1665|74|0.184539Å|
| Hans8FNR | A24–133,110residues |1665|179|0.174741Å|
| Mex8FNS | A29–133,105residues |1544|167|0.229532Å|

Original C–H/N–H/O–H medians were around1.18–1.20Å, like the old1H4I problem.
The existing whole-chain preparation projects each protein H along its actual
parent direction to its ff19SB length; all audited residuals now fall below
7e−15Å. Water O–H lengths are0.9572Å. Heavy atoms have zero physical displacement.
These are radial repairs, **not optimized hydrogen orientations or a proton
ensemble**. Existing water H–O–H angles retain their original103.30–106.90° range.

All three protonation receipts declare pH5, no missing-residue reconstruction,
no missing-heavy additions and no terminal completion. Complete actual peptide
graphs have109/109/104 bonds. The8FNR missing-residue annotation concerns a
nonselected chain. No nonbonded heavy pair below1.8Å remains in the selected
protein/water sources after excluding1–2 and1–3 graph neighbors. This narrow
geometric screen does not establish equilibrium or rule out milder clashes.

Source donor bonds are intact:

| Source | Asp85 CB–CG (MexAsp86) | Glu91 CG–CD (MexGlu95) |
|---|---:|---:|
| Hans8DQ2 |1.516147Å|1.520734Å|
| Hans8FNR |1.534019Å|1.540154Å|
| Mex8FNS |1.513204Å|1.529219Å|

The >2.7Å covalent failures belonged to rejected whole-chain proposals, not these
origins. MexAsp92 CB–CG is1.524536Å; Mex residue92 is Asp, not HansGlu91.
Do not mechanically equate residue numbers across proteins.

## Exact proposed local region

Use the same deterministic rule in each source before any new scores:

1. Entire EF3 loop: HansA83–94 or MexA84–95, including every real heavy atom and
   normalized H. This retains complete donor sidechains, local peptide carbonyls,
   nearby charged groups and the internal hydrogen-bond scaffold.
2. Preserve the bounding peptide units by adding preceding C/O (Hans82,Mex83)
   and following N/H (Hans95,Mex96). Cut only preceding CA–C and following N–CA;
   source-mapped caps are not independent physical degrees of freedom.
3. Retain the nonlocal support amide: HansA65 C/O plus actual bonded A66 N/H;
   MexA66 C/O plus actual bonded A67 N/H. Their carbonyl O lies2.752/2.799/2.832Å
   from the core backbone N90/N91. Keep the actual C–N bond; cap only C–CA and
   N–CA sigma attachments. These atoms are selected by actual connectivity,
   not merely residue-number adjacency.
4. Include all original direct-core waters and any actual water O within3.4Å
   of the complete original-core N/O inventory. Use complete water O/H/H.
   The selection is a polar-contact rule, not proof every water forms a hydrogen
   bond in its inherited orientation. Freeze membership during perturbation.

| Source | QM atoms including exchanged metal and4caps | Selected source water residues | Proposed charge |
|---|---:|---|---:|
| Hans8DQ2 |195|A337,A350|−1|
| Hans8FNR |198|A314,A345,A369|−1|
| Mex8FNS |209|A320,A331,A358,A373,A376,A394,A434|0|

The charges follow the actual current residue templates: Hans loop five acidic
sidechains plus Lys84 gives−4; Mex loop five acidic sidechains plus Lys93/Lys94
is−3. Neutral added amide/water groups andLn(III) give the stated proposals.
Preparation must still verify exact electrons/parity/cap valence. La physical
multiplicity1 and Dy6 are state hypotheses requiring a genuinely Dy-capable
reference; native f-in-core singlet support is not that validation.

`RESULT_v2.json` lists every selected physical index/source ID, all four boundary
attachments, full water inventory and source pins. `RESULT.json` is the preserved
initial contact inventory before the nonlocal support amide was closed; use v2.
The old graph is reusable for chemical identity, but old XYZ H coordinates differ.
These larger regions intentionally qualify local chemistry/support, not screening
throughput. Measure one reference before expanding the finite matrix.

## Occupancy: retain what is actually deposited

The first local fixed-environment experiment should retain every actual source
spectator, identically for targetLa and targetDy:

| Source | EF1 | EF2 | EF4 | EF2–EF3 distance |
|---|---|---|---|---:|
| Hans8DQ2 |La(+3)|La(+3)|**Na(+1)**|11.969Å|
| Hans8FNR |Dy(+3)|Dy(+3)|Dy(+3)|11.831Å|
| Mex8FNS |Nd(+3)|Nd(+3)|Nd(+3)|11.410Å|

EF1/EF4 lie about29.9–30.4Å from EF3. Use explicitly declared frozen formal
monopoles for this electronic scout only, not invented MM atom types/LJ/spins.
This approximates spectator electrostatics without their polarization or local
quantum response; it is not a relaxation model. DQ2Na is not replaced byLa/Dy.
Do not reuse EF23/EF1234 endpoint XYZ files blindly: those are separately
constructed occupancy hypotheses, some deleting spectators or replacingNa.
The source `atoms.json` excludes metals by design; append every actual spectator
from the pinned source site inventory and place only EF3 in QM.

With protein charges−4/−4/−10 and neutral water, the actual all-site systems
have charges+6/+8/+2. Given the proposed QM charges, consistent MM/field totals
would be+7/+9/+2. This is a closure target for explicit local boundary accounting,
not permission for global neutralization or ad hoc charge redistribution.

## Reusable export and missing implementation

`scripts/affordable_state.py:source_protein` is not directly usable unchanged:
its wet-chain input check rejects HOH, it has no spectator handling, and it reads
coordinates from the old protonated PDB. It does not itself add/normalize H, so
using its returned coordinates would reintroduce the documented defect.

Existing whole-chain audits retain topology and total template charge but no
exported per-atom embedding charge array (`ff_atom_charges_used_for_scoring=false`).
A small new research export is therefore needed:

- Build ff19SB+TIP3P topology/System from actual selected-chain protein/water,
  explicitly excluding only separately recorded spectator/target ions from FF
  typing. Do not call addHydrogens, repair missing atoms, or optimize.
- Map charges to **normalized** `atoms.json` coordinates by the exact recorded
  identities; independently confirm atom inventory and charge sums.
- Define local cap/charge-shift accounting for the four CA–C/N–CA cuts. Preserve
  real MM boundary atoms in the physical mapping. Fractional classical charges
  removed with a capped QM fragment need not equal its integer formal charge;
  show the local redistribution and charge closure explicitly. The old1H4I
  sidechain-only redistribution should not be copied unchanged across these
  peptide boundaries.
- Add the declared spectator monopoles without double counting them in FF or QM.
  Keep gauge, field, occupancy, proton inventory and all external coordinates
  identical across targetLa/Dy. Retain native charge and model provenance.

Thus the preparation is feasible but no ready-to-run local field presently exists.
No molecular parameterization or charge export was executed during this audit.

## One bounded local perturbation proposal

Use +2° right-hand rotation around actual CA→CB of HansAsp85/MexAsp86. Rotate
the **whole distal sidechain**, including CG,OD1,OD2 **and HB2/HB3**, keeping
all other atoms fixed. Rotating only the carboxylate heavy atoms would distort
CB–H versus CB–CG angles. Source full-graph mapping and contact checks must pass
before admission; propose maximum physical displacement0.15Å, no donor deletion,
water change or proton change. The support-water inventory and spectator field
remain fixed. Apply the exact same A/B pair to both target metals. This probes
nearby source-derived restoring response, not preference, population or whole
protein accommodation. Glu91 remains complete in the electronic region and its
force projection is retained as a diagnostic.

Root owns the declared reference Hamiltonian, electronic-state qualification,
actual preparation/execution and any later expanded-boundary test. No old failed
LanM campaign, within-series threshold or reserved-library outcome is reused.
