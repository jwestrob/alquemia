# Asn290 discrete-conformer preparation review

Status: preparation and conditional scoring implementations passed bounded static
review. Actual Slurm preparation/admission remains required.
This review uses existing code and saved metadata only. No scientific imports,
geometry calculations, tests, energy calls or submissions were performed.

## Material contract

The authority is `../ASN_ROTAMER_PLAN.md`, using the completed inventory
`../a0ac_asn_inventory/EVIDENCE.json`. The inventory records the exact204-atom
state and source graph for all ten A0ACD6B9F2 origins, with neutral Asn290 and
both physical HD21/HD22 hydrogens attached to ND2. Caps terminate outside the
moving terminal group. No new inventory or geometry survey was run here.

A proper180° rotation of OD1, ND2, HD21 and HD22 about each actual CB→CG axis
is a chemically defined conformer. CG stays on the axis; atom identities are
retained. This is not O/N relabeling and does not establish the experimentally
correct rotamer. The existing correctly called La-conditioned sample4 remains
an essential counterexample in the complete ten-source denominator.

Preparation must validate each named source/QM map, exact graph and H parents,
rigid-amide and affected-bond invariance, fixed CA/backbone stereochemistry,
fixed caps/anchors and paired Ca/La coordinates. The declared discrete flip is
outside the old0.8rad/0.8Å donor-search domain; no old admission result may be
silently relaxed or passed retrospectively.

## Full-source contacts

Use all saved physical source heavy atoms and hydrogens with the inherited
relative severe-contact predicate: new distance<1.0Å heavy/heavy or<0.55Å with
H, while the same unrotated pair was at or above its limit. Retain inherited
short contacts separately. Passing this narrow gate does not prove a continuous
collision-free rotation or a relaxed equilibrium structure.

The four moving source atoms must appear exactly once, rather than remain as
stale duplicate origin centers in the exterior. Curie's agreed construction
uses the local204 block followed by only source identities absent from its real
atom map. Three unchanged modeled PQQ hydrogens remain explicit; synthetic caps
are distinguished from physical source atoms. Any native La in a La-conditioned
source requires explicit parser support, without relabeling its actual source
identity or guessing an unmapped atom.

Failures must retain the ten-source/twenty-endpoint denominator. No change of
angle, proton inventory or source selection may rescue an inadmissible flip.

## Final preparation implementation

No remaining concrete static blocker was found in
`../a0ac_asn_rotamer/PREPARE.py` SHA256
`45edb9a955782b67b07f8627e42e9e05a4f92c12fb2112fa1a620bcb053e111c`
and `RUN.sbatch` SHA256
`e7863329376d137552765b73a58e4ea332edf73d7d94ecb9d79e208db3558273`.
The reviewed plan SHA256 is
`c053d49f063bd3bc383766c1223e0996f784c72539fec005ca38bd4d915ee1f1`.

The source verifies the actual CB–CG/CG–OD1/CG–ND2/ND2–H graph and exact source
hydrogen parents. It independently constructs a proper180° rotation and compares
it with the existing physical kinematic map, without invoking the old local
search admission domain. All four amide atom identities remain intact. It checks
rigid and covalent invariants, fixed stereocentre/substituents and fixed caps.

The contact union deduplicates by source identity and includes all physical
source hydrogens. Separate nearest omitted-heavy and H-involving records retain
before/after distances and identities; inherited short contacts and synthetic-cap
contacts remain distinguishable. The parser explicitly permits native source La
without renaming atoms. Unsupported rows retain their source identity and fail
closed. Actual numerical preparation is still required before energy dispatch.

One auxiliary replay issue was corrected during review: kinematic rotation can
introduce roundoff at CG even though CG is on the axis. The final code requires
the written CG coordinates to stay exactly unchanged and bounds only this
auxiliary replay residual by the already declared1e−10Å invariant tolerance.
This does not loosen the old donor domain or move the emitted axis point.

## Conditional scorer review

At Gauss's request, a bounded static read also covered
`../a0ac_asn_scoring/SCORE.py` SHA256
`7682a6265084039a1acc663bf4812cfff3aebd7b3acbf454a0f6642ab1701388`
and `RUN.sbatch` SHA256
`3bc0ccf32d6424fde3f1fa7fb5d33d15a5137c64c942ad58b9865cf0bf040aa1`.
No material static blocker was found. Its expected source and paired-cell schema
matches the final preparation implementation. Actual geometry output still needs
to pass and be joined before execution.

The shared scientific predicate explicitly checks204 atoms, the exact state,
model, coordinates and force convention. Archived CPU and H200 receipt schemas
are handled separately; H200 execution is not relabeled as seven-thread CPU.
All sixty real archived force arrays are checked for finite204×3 shape inside
Slurm. New cells retain strict seven-thread CPU, version, unchanged-parameter
and single-forward admission from the qualified worker.

Every preparation source remains visible, including unsupported flips. Both
metals must be admitted before a source is dispatched. Complete shared common4
selection keeps the original three first in tie order and reports mathematical
and strict−0.1kcal/mol operational choices separately. A missing new cell leaves
common4 unavailable; all original common3 baselines remain separately named.
Matched-coverage spreads identify their actual source set. The sign-blind
absolute1kcal/mol continuation gate uses only complete source comparisons and
does not launch solvent automatically or assign new labels.

A concrete path-only provenance mismatch was corrected before freeze: preparation
binds the workspace inventory, while the initial scorer expected the diagnostic
copy's entire pin dictionary. The scorer now retains the actual preparation pin
and verifies both copies against the identical fixed content SHA. It also
preserves partial results and execution receipts through cleanup/reporting
failures. Neither change alters inputs, physics, selection or admission.

## Actual geometry review — job4614

The saved preparation admits all ten sources and twenty exchanged-metal cells.
No omitted-source packing obstruction is evident in the recorded nearest contacts.
This supports the declared fixed-conformer energy test; it does not establish a
relaxed alternate rotamer, a collision-free rotation path or a complete protein
energy. The unchanged severe-contact gates remain the admission criteria.

Evidence: `workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_preparation_v1/execution/PREPARATION.json`,
SHA256 `508e7594aea7ad6198c62ce931d84a5f19d72b71070edcdbc782787f8cd9ca30`,
and each row's `GEOMETRY_ADMISSION.json` and `CONTACT_IDENTITIES.json`.
The manifest SHA256 is
`ae4184400f6908039ebb1eaa2121181e92510738d88366d5175b9ecc9df624fb`.
Stored nearest-contact identities agree with their saved union-index identities.
This review reads saved results only; it performs no distance or energy calculation.

Each receipt includes9529 full-source atoms and9342 appended physical atoms,
with all physical hydrogens included and synthetic caps distinguished. All ten
record no new severe contacts and no inherited contacts below the frozen severe
thresholds. Covalent, rigid-group and auxiliary replay residuals are about1e−15Å;
the emitted CG, other coordinate lines and cap anchors remain fixed. The roughly
2.3Å terminal-group displacements belong to this explicit discrete conformer
protocol, not the old0.8Å local search domain.

Nearest omitted-source distances below are copied from the saved receipts and
rounded for display. Ca/La here denotes source folding conditioning, not a new
biological label. All listed atoms are on source chainA.

| Source | Nearest omitted heavy contact | Before→after (Å) | Nearest omitted H contact | Before→after (Å) |
|---|---|---:|---|---:|
| Ca0 | Asn290 ND2—Gly289 CA | 6.394314→4.505337 | Asn290 HD21—Gly289 HA2 | 6.735458→3.052105 |
| Ca1 | Asn290 ND2—Gly289 CA | 6.440944→4.572510 | Asn290 HD22—Gly289 HA3 | 6.704466→3.464226 |
| Ca2 | Asn290 ND2—Gly289 CA | 6.348390→4.475576 | Asn290 HD21—Gly289 HA2 | 6.822282→2.741161 |
| Ca3 | Asn290 ND2—Gly289 CA | 6.407013→4.519666 | Asn290 HD21—Gly289 HA2 | 6.849887→2.867103 |
| Ca4 | Asn290 OD1—Gly289 CA | 6.147383→5.232481 | Asn290 HD21—Gly268 HA3 | 6.029638→3.604475 |
| La0 | Asn290 ND2—Gly289 CA | 6.391011→4.518698 | Asn290 HD21—Gly289 HA2 | 6.731266→3.034289 |
| La1 | Asn290 ND2—Gly289 CA | 6.435257→4.590598 | Asn290 HD22—Gly289 HA3 | 6.794588→3.356908 |
| La2 | Asn290 ND2—Gly289 CA | 6.353987→4.512834 | Asn290 HD21—Gly289 HA2 | 6.843647→2.817094 |
| La3 | Asn290 ND2—Gly289 CA | 6.400332→4.551274 | Asn290 HD21—Gly289 HA2 | 6.821932→2.917754 |
| La4 | Asn290 OD1—Gly289 CA | 6.106476→5.274727 | Asn290 HD21—Gly268 HA3 | 5.776476→3.899380 |

One distinct local limitation matters for interpreting energies: several flips
place Asn290 HD21 close to its backbone H, both already inside the204-atom scored
core. Examples are Ca2 at0.961077Å and La2 at0.960025Å, compared with their original
4.701552Å and4.680432Å separations. These are nonbonded H–H contacts, unlike the
nearest H records in sample1/4 that identify the preserved covalent ND2–H bond.
They pass the existing0.55Å severe-contact gate; no stricter threshold is introduced.
Their energetic cost is included locally. An unfavorable flip energy can therefore
reflect this unrelaxed conformer's strain and does not rule out every relaxed
alternate rotamer basin. The omitted-contact records do not support attributing
such a penalty to an unmodeled exterior clash.

The actual one-CPU preparation receipt reports5.166875s wall,4.195694 CPU-s,
93748KiB peak RSS and zero energy calls. No additional preparation, optimization,
submission or scoring modification was performed for this review.
