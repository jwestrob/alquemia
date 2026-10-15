# Intact-PQQ pose: bounded static review

Status: mapping and inherited-gate review complete; implementation review pending.
This review uses source code and saved metadata only. No scientific imports,
coordinate calculations, tests, model calls or submissions were performed.

## Coordinate and chemistry contract

The current authority is `../PQQ_POSE_PLAN.md`: two sources, three saved donor
bases, both signs of a one-degree rotation and both metals (24 new cells).
The twelve original cells remain the baseline; the complete common pool has
nine geometries per source. The earlier sixteen-cell suggestion in
`RESULT_REPORT.md` is superseded by this explicitly broader plan.

Use the source's named PQQ map, including all 24 heavy atoms and three physical
hydrogens. In the generic schema, e1 follows C5→O1, e2 is the perpendicular
component of C5→C6, and e3=e1×e2. The center is the 24-heavy centroid.
H_N1 is attached to N1, H_C3 to generic C10, and H_C8 to generic C1.
The generic names must not be silently replaced with CCD numbering.

The inherited donor map includes PQQ heavy atoms but omits the three modeled
PQQ hydrogens. Replaying original donor admission and then explicitly checking
the intact 27-atom rigid transformation is appropriate. This new motion must
not be represented as an old donor-mode replay. All internal distances and
physical H-parent distances, fixed atoms, caps and cap anchors require the
declared checks.

## Actual inherited overlap rule and full-source reuse

`scripts/accommodation_torsion_profiles.py::geometry_check` rejects a pair only
when its new distance is below 1.0 Å (two heavy atoms), or 0.55 Å (either atom H),
and its original distance was at or above that same limit. It does not invent
new bonded exclusions. Here the comparison baseline is the corresponding saved,
unrotated donor candidate, as the pose plan declares.

The existing source unions allow the same check against omitted physical atoms:

- `workspaces/native_cavity_shielding_20261014/preparation_v1/execution/PREPARATION.json`
  SHA256 `501105d36b69ea915e93be4a082c9bda3b40bdb61275019cb05226e59cf9e460`:
  four prepared origin unions.
- `workspaces/native_cavity_accommodation_20261014/preparation_v1/execution/PREPARATION.json`
  SHA256 `3b50b7514e8b306af11deb0a33f73dd61513a860fdaec81afbd40e92315aa97f`:
  eight prepared donor-candidate unions.

Each contains 9,474 centers, with the exact local 192 first and unchanged source
appendages thereafter. Pin the corresponding source/state/map/XYZ, verify that
local block against its saved candidate, and test rotated PQQ against appended
atoms. The existing accommodation preparation provides this relative-contact
implementation. Synthetic caps are not physical source atoms. No environment
energy, reconstructed source geometry or new union preparation is needed.
Passing this extreme-overlap check establishes no new severe collision; it does
not establish general steric relaxation or mechanical support by the protein.
