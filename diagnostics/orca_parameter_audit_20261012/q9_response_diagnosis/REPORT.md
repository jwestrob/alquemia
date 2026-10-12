# q9 saved-motion diagnosis: local Glu response remains unresolved

12 October 2026. Coordinate and archive inspection only; zero molecular calls,
optimizations, submissions or production changes. The reported ~12.1 kcal/mol
Ca−La accommodation shift motivated this audit; this directory does not recompute
that energy or assign its physical cause.

**The dominant common motion balances the two Glu213 metal–oxygen distances,
tightens its existing Trp304 H-bond geometry, and loses its existing backbone
H-bond geometry. Coordination number and denticity remain unchanged. Both saved
La proposals reach the 0.8 Å displacement boundary.** Neither the original nor
expanded context includes nearby Leu485 packing, but the omitted atoms produce
no new heavy pair within the predeclared 3.5 Å descriptive bin. Geometry does not
establish missing steric resistance as the cause of the energy response.

## Exact scope and definitions

The source is `q9z4j7-pqq-la_model__conditioned_La__seed-1_sample-0` from the pinned
`context_factorial_preparation_v1/PREPARATION.json`. It supplies the actual origin,
`original_adaptive_La` and `expanded_adaptive_La` saved proposals, both context
maps and proposal receipts. Original/expanded models have 145/192 atoms and
9/14 artificial caps; charges are Ca−3/La−2, singlets. Source-mapped motion is
identical when each saved proposal is grafted into either context. Metal, PQQ
and other real scaffold coordinates remain fixed. Thirteen real atoms move,
seven heavy, all in Glu213/Asp350.

Definitions were recorded in [PLAN.md](PLAN.md) before coordinate inspection:
existing chemically typed metal donors ≤3.1 Å; diagnostic donors ≤4.0 Å;
descriptive contact bins ≤3.5, 3.5–4.0 and >4.0 Å; N/O–H···N/O geometric flag
requires donor–acceptor ≤3.5 Å, H–acceptor ≤2.5 Å and angle ≥120°. These bins
are not fitted energy, clash or affinity thresholds. All contacting real source
atoms are retained, regardless of whether they are favorable or unfavorable.

The full-parent coordinate graft is diagnostic and has not been energetically
evaluated. SourceGraph/ff19SB bonds identify covalent neighbors; pairs separated
by one/two bonds are excluded from nonbonded interpretation, three-bond pairs
are explicitly marked as 1–4. All source bonds involving moved atoms are
preserved to 1.6×10⁻¹⁵ Å. Artificial caps never count as real protein contacts.
`RESULT_V2.json` supersedes the first result by including same-residue H-bonds
and source-bond exclusions; all distance/angle thresholds are unchanged.

## Metal coordination and actual torsions

All three geometries have **CN7: six O plus one N**. Glu213 is bidentate,
Asn300/Asp350 monodentate, PQQ tridentate. There is no donor gain/loss or
denticity switch under the existing classifier.

| Donor | Origin Å | Original La proposal Å | Expanded La proposal Å |
|---|---:|---:|---:|
| Glu213 OE1 | 2.254 | 2.680 | 2.656 |
| Glu213 OE2 | 2.940 | 2.600 | 2.620 |
| Asn300 OD1 | 2.370 | 2.370 | 2.370 |
| Asp350 OD1 | 2.929 | 2.876 | 3.035 |
| Asp350 OD2, diagnostic only | 3.708 | 3.405 | 3.626 |
| PQQ O1 | 2.361 | 2.361 | 2.361 |
| PQQ N2 | 2.438 | 2.438 | 2.438 |
| PQQ O6 | 2.595 | 2.595 | 2.595 |

Actual changes computed from four real source coordinates agree with the saved
rotation parameters. Their signs use the source-map atom-order convention.

| Torsion | Original La change, degrees | Expanded La change, degrees |
|---|---:|---:|
| Glu213 χ1 | −10.121 | −10.129 |
| Glu213 χ2 | +3.525 | +4.266 |
| Glu213 χ3 | +20.324 | +19.885 |
| Asp350 χ1 | +10.945 | 0 |
| Asp350 χ2 | 0 | −9.663 |
| Asn300 χ1/χ2 | 0/0 | 0/0 |

Glu213 OE2 is the boundary-setting atom: displacement 0.800000000275 Å in the
original proposal and 0.800000000010 Å in the expanded proposal. Both receipts
say `boundary_flag=true`, one boundary heavy atom, no angular boundary, physical
mapping feasible, no new severe overlaps, and no unconstrained-minimum claim.
Their saved vacuum-OMOL omitted angular loads remain 10.145 and 12.136
kcal/mol/Å, respectively. The selected-mode maxima are 1.895/1.692. These are
constrained proposals, not demonstrated stationary physical minima.

## Real second-shell changes

The strongest common geometric changes involve groups **already included in
both contexts**:

| Donor–H···acceptor | Origin D–A / H–A Å | Original La D–A / H–A Å | Expanded La D–A / H–A Å | Geometric flag |
|---|---:|---:|---:|---|
| Glu213 N–H···OE2 | 3.050 / 2.193 | 3.572 / 2.874 | 3.582 / 2.880 | Lost in both |
| Trp304 NE1–HE1···Glu213 OE1 | 3.156 / 2.298 | 2.784 / 1.865 | 2.804 / 1.891 | Retained, tighter |
| Arg378 NH2–HH22···Asp350 OD1 | 3.438 / 2.552 | 3.895 / 3.008 | 3.291 / 2.382 | New in expanded proposal only |

Glu backbone angles change 126.5→117.0/117.3°; Trp angles 127.1→130.7/130.3°.
The Arg origin narrowly fails the fixed H-distance threshold; its expanded
angle is 132.3°. These are geometry flags, not hydrogen-bond energies. Ser352
backbone N/H is expanded-only: its D–A distances 3.750→3.954/3.615 Å fail the
fixed 3.5 Å criterion throughout. Asn300 ND2 heavy contacts do not meet the
explicit-H angle criterion and are not counted as H-bonds.

At the same 3.5 Å heavy-pair boundary, original La creates only Asp350 CG–Trp304
NE1 (3.640→3.428 Å); it loses Glu213 OE2–PQQ C13 (3.096→3.556), Asp350
CG–Arg378 NH2 (3.279→3.510), Asp350 OD1–Arg378 NH2 (3.438→3.895) and Asp350
OD1–PQQ O1 (3.432→3.686). Expanded La creates Asp350 OD1–Ser352 CB
(3.526→3.444), and loses Glu213 OE2–PQQ C13 (3.096→3.562) and Asp350
OD1–Asn300 ND2 (3.499→3.658). All partners in these changes are already in
both contexts. Carbon-containing contacts are not automatically H-bonds.

## What both contexts omit

There are 28 omitted source atoms within 4 Å of at least one moved atom in at
least one geometry: 13 heavy and 15 H. This count includes ordinary backbone
proximity and weak peripheral H contacts; it is not a list of 28 missing forces.
All 51 moved–omitted pairs are in [OMITTED_CONTACTS.csv](OMITTED_CONTACTS.csv),
with identical distance bins, actual source-bond separations and 1–4 flags.

| Omitted interaction | Origin Å | Original La Å | Expanded La Å | Interpretation |
|---|---:|---:|---:|---|
| Glu213 CG–Leu485 CD2 | 3.836 | 3.599 | 3.599 | Largest common omitted heavy compression; stays above 3.5 |
| Glu213 HG3–Leu485 CD2 | 2.873 | 2.585 | 2.594 | H–C packing tightens; no energy evaluated |
| Glu213 HG3–Leu485 HD21 | 2.717 | 2.415 | 2.427 | H–H packing tightens |
| Asp350 HB2–Phe442 HB3 | 2.488 | 2.281 | 2.488 | Original-only H–H compression |
| Asp350 OD2–Phe442 CD2 | 3.714 | 4.015 | 3.773 | Moves away in both proposals |
| Asp350 OD2–Phe442 CB | 3.893 | 4.204 | 3.991 | Moves away in both proposals |
| Asp350 HB2–Trp349 O | 2.536 | 2.440 | 2.536 | Backbone C–H···O proximity, not N/O–H flag |

No omitted heavy pair is ≤3.5 Å at origin or either endpoint. At that boundary
the omitted-pair inventory has **zero newly formed pairs**, heavy or H, and one
lost H–H pair, Glu213 HB3–Trp591 HH2 (3.465→3.557 Å in both). At 4 Å, Glu213
HG2–Leu485 HD23 and HG3–its own backbone O newly enter the 3.5–4 Å bin; several
peripheral pairs leave it. No omitted donor/acceptor loses or gains a qualifying
N/O–H H-bond. The short Glu HB–own C and Asp HB–own N contacts are marked 1–4,
not presented as new clashes. Leu485 packing is a concrete missing-environment
candidate; the geometry gives no quantified resistance or unique causal result.

## Why expansion still excludes these atoms

The existing builder selects **N/O/S neighbors of fixed donor anchors**, then
complete sidechain/peptide fragments with sigma caps. Original discovery radius
is 3.5 Å; expanded radius is 4.3 Å. It is neither an all-heavy packing shell nor
a swept-motion selection. The q9 expanded fragment union came from La-conditioned
samples 1–4 and was transferred to canonical sample0. The replay in
[OMISSION_POLICY.json](OMISSION_POLICY.json) verifies every relevant recorded
contact against those four actual source structures.

- Leu485, Ala483, Phe214 and Phe442 sidechains have no N/O/S selector. Their
  carbon/H packing does not independently trigger selection.
- Glu213 O and Phe214 N, which would select their intervening peptide, remain
  5.735–5.796 and 5.325–5.421 Å from the nearest fixed anchor across the four
  source members. Selected Glu sidechain/preceding peptide does not include all
  of Glu's backbone.
- Trp349 O is 4.889–4.958 Å away; Asp350 N is narrowly outside at 4.310–4.324 Å.
  Their peptide is omitted. Asp350's selected sidechain/following peptide does
  not include the preceding peptide.
- Trp591's polar NE1 selector is 7.258–7.364 Å away although some ring/H atoms
  lie within 4 Å of moved Glu hydrogens.
- Phe442 O **is selected** at 3.969–4.084 Å in expanded discovery. Its peptide
  terminates through the existing capped fragment boundary; this does not bring
  in its nonpolar sidechain or CA/HA.

These facts explain the actual omissions; they do not justify adding all 28
nearby atoms or assert that a radius change will fix the energy response.

## Archived DFT compatibility

[DFT_INVENTORY.json](DFT_INVENTORY.json) records exact source pins and raw-output
checks for the completed q9 Ca/La origin pair: r2SCAN-3c/CPCM(Water), DefGrid3,
73 atoms, charges **Ca−2/La−1**, singlets. Input/output/XYZ/manifest/receipt hashes
match, both outputs terminate normally, and each archived XYZ exactly matches
the current frozen-core origin coordinates. These are actual completed endpoints.

They are **not compatible origin-to-proposal work measurements** for the present
145/192-atom, Ca−3/La−2 models. The bounded related inventory contains only the
completed 4MAE/two-PLM extra-Asp ±0.2 rad DFT paths; that experiment explicitly
did not test Glu/coupled or relaxed candidates. The later strict-response study
reuses those 18 points. The separate native-candidate campaign remained dry-run
only after its gate failed. No exact q9 motion DFT pair was identified in these
named sources; this is not an exhaustive claim about every historical archive.

## Interpretation and retained evidence

The common response can be described as Glu carboxylate balancing plus a change
between two already-modeled H-bond geometries. The expanded proposal also changes
the active Asp torsion and gains an Arg H-bond geometry, so it is a distinct pose.
Neither a CN increase nor an unmodeled short heavy clash explains the result.
Modest omitted nonpolar compression and local model energetics remain separate
hypotheses. The saved 0.8 Å bound and residual angular loads also limit any
minimum/accommodation interpretation. Geometry and the incompatible archived
DFT origins cannot distinguish these hypotheses; no new DFT is proposed or run
by this report.

Primary evidence: [RESULT_V2.json](RESULT_V2.json), its pinned source maps/XYZs,
proposal receipts and ff19SB source graph. Analysis implementation and source
hashes are recorded in the JSON files. `ANALYZE.py`/`RESULT.json` are retained
as the earlier pass; `ANALYZE_V2.py`/`RESULT_V2.json` are the report's authority.
`SUPPLEMENT.py` independently replays the four source-member selector distances,
verifies the finite DFT artifact set and writes the contact CSV. No archive-wide
hash sweep, model invocation or numerical energy evaluation occurred.
