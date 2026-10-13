# Recovery of fifteen exact-PQQ preparation failures

Declared before new coordinate inspection or scores, 12 October 2026.

Scope: all five saved Ca-conditioned seed-1 samples of C5AXV8, Q88JH5 and
Q9Z4J7. Preserve the original 125-source result and frozen canonical bands.
The old C5AXV8 audit already finds intact PQQ with Ca about 25 A away; an
off-site ion is not repaired by renaming or moving atoms.

The finite job will inspect the exact pinned raw, normalized and protonated
structures for these fifteen sources. It will replay the original schema-aware
normalizer on the saved residue-alias CIF. It will report PQQ identity/formula,
the exact selected metal, nearest PQQ and protein O/N distances, the existing
typed donor contacts, and all-heavy-atom identity/coordinate preservation.
Distance gates remain those in the original code: PQQ association 4.0 A,
typed coordination 3.1 A, at least six compatible typed donors for the fixed
core. No score, experimental class or desired prediction selects a geometry.

Only a demonstrated metal/PQQ residue-name discrepancy between the saved
protonated file and the unchanged normalizer replay can receive a metadata-only
repair. Chain, residue number, insertion code, atom name, element, altloc,
occupancy and every coordinate remain unchanged. Existing hydrogens are reused;
there is no reprotonation, new water, cap adjustment, atom relocation or refold.
A chemically valid but distant PQQ gives `requires_new_site_conditioned_structure`.
Missing/invalid/ambiguous PQQ gives an explicit unavailable reason.

Any associated source proceeds through the unchanged fixed-core builder and
complete context builder, followed by its already frozen 4.3 A canonical
fragment union and original state-anchor checks. No chemical gate is bypassed.
Prepared rows/maps/endpoints use the existing preparation schema for later root
scoring. No MACE, xTB, ORCA, optimization or molecular dynamics calls occur.
All coordinate parsing, preparation and scientific checks run in Slurm with one
numerical thread. Login work is limited to code and metadata/byte pinning.

Outputs: immutable per-source AUDIT.json and ADMISSION.json, PREPARATION.json,
AUDIT.json, REPORT.md and EXECUTION.json under the manifest's fresh output path.
All fifteen remain in the output even if no repair is justified. The next step
for an off-site source is a separately declared structure-regeneration protocol;
this job does not fabricate an on-site structure.
