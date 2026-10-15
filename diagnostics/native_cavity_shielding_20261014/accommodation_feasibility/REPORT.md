# Accommodation shielding is feasible from eight exact saved cells

**All eight nonorigin cells are available; none requires new SCC, MACE, DFT,
optimization or folding.** The same-source union can retain its identity/order
while replacing the first 192 coordinates and charges with each saved candidate.
This still needs a Slurm preparation/admission check; no accommodated union was
constructed or scored by this audit.

## Exact reusable evidence

[INVENTORY.json](INVENTORY.json) contains all 12 local cells (four origins plus
eight nonorigins), exact XYZ/maps/state/source pins, raw vacuum charge and
execution receipts, native component results, original component-v2 acceptance,
and the completed origin-union artifacts. SHA256:
`9b38d4b184ddb69079fbd955bc517ec4abb642d1d71b011b734cd7efb01d7658`.
The metadata-only [audit](INVENTORY.py) verified 528 artifact hashes; it performed
no coordinate arithmetic or energy analysis.

| Source | Saved donor arrangement | Ca | La |
|---|---|---|---|
| ExaA two-Ca sample0, control | adaptive_Ca | Available | Available |
| sample0 | adaptive_La | Available | Available |
| ExaA two-Ca sample3, target | adaptive_Ca | Available | Available |
| sample3 | adaptive_La | Available | Available |

All retain the same 192 atoms, 14 caps, dry oxidized-PQQ3− hypothesis, Ca−3/La−2
singlets and 588 GFN valence electrons. Each candidate has its **own exact
vacuum populations**; origin charges must not be substituted. Every native cell
has passed water-diagonal/precision checks and a separately accepted
standalone-xTB6.7.1/tblite0.7.0 solvent-component pair. The historical failed
absolute crossbackend gate remains failed.

Authoritative source and result roots:

- `workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json`
  (SHA `757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343`).
- `workspaces/ca_coverage_recovery_20261012/two_ca_scoring_q9_v1/`:
  actual proposal coordinates, cross-metal cells, proposal manifests and maps.
- `workspaces/ca_coverage_recovery_20261012/two_ca_solvent_q9_v1/`:
  qualified per-geometry vacuum/water outputs and component collection.
- `workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/transfer_v1/tasks/`:
  exact local fixed-vacuum-charge native components.
- `workspaces/native_cavity_shielding_20261014/preparation_v1/execution/PREPARATION.json`:
  each origin's 9,474-center capped union and 9,460-center cap-removed map;
  `native_v1/execution/FINAL_COLLECTION.json` contains their completed energies.

## Motion and boundary feasibility

All four saved proposals select the same physical modes: Glu213 χ1/χ2/χ3 and
Asp350 χ1 (ordering differs). Their mode membership contains 15 physical atoms,
all mapped into the local block: Glu213 CB/CD/CG, HB2/HB3/HG2/HG3/OE1/OE2;
Asp350 CB/CG, HB2/HB3/OD1/OD2. Axis endpoints may remain stationary despite
membership. No selected moving member is absent from the local block, and
**none of the 14 cap links touches a selected moving member**. The fixed CA
axis/flanking atoms are present in the original physical map. Metal translations
are zero. Existing proposals passed their source/cap/bond/overlap checks;
sample3 adaptive_La retains its recorded displacement-bound flag.

Thus appended physical atoms can remain at their exact original coordinates;
no missing flanking atom has to move, and no new cap construction is needed.
However, the old physical map covers 189 atoms, not the entire protein. It does
not qualify new close contacts between the moving donors and the newly included
source atoms. Existing saved analyses identify Leu485/Phe442 packing contacts;
they do not replace the required full-union admission.

Before any native call, replay each saved q against its pinned original
kinematics, match the candidate XYZ at the existing tolerance, and verify:

- Paired Ca/La geometry, atom/state order and exact candidate charge provenance.
- Only mapped local coordinates change; appended source coordinates and exact
  zero charges remain fixed, with no duplicate original-position copies.
- Original physical displacement/bond/cap gates and all fixed atoms are retained.
  Check new severe overlaps against the origin union using the existing
  1.0 Å heavy / 0.55 Å H-involving policy; separately preserve known synthetic
  cap/source contacts rather than mislabeling them as intact chemical geometry.
- Cap-removed controls remove the same 14 indices and retain the same real
  candidate charges. No radius edits or charge redistribution occur.

This is a finite test of **accommodation-dependent electrostatic exclusion**,
not a revival of the unsuccessful origin-only remedy. Retain the complete
origin/adaptive_Ca/adaptive_La pool for both sources/metals and all failures.
Keep native HB, SASA and reference terms separate; no calibrated cross-size
classification or source-only physical-cavity claim follows. The completed
origin result's aggregate precision flag remains a limitation, not permission
to loosen the new contrast's numerical gates.
