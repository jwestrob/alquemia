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
