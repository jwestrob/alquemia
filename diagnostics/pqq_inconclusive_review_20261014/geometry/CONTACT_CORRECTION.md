# Saved contact descriptor correction

Job4038 completed all48 cases. Its new nearest-nonbonded descriptor was invalid:
the kinematic bond graph omits fixed PQQ internal bonds, so real cofactor C–O
bonds appeared in the ranked list. Preserve that output and do not interpret its
`closest_mapped_nonbonded_heavy_pairs` or `nonbonded_scope` fields. The donor,
confidence, PQQ, motion and recorded optimizer fields are unaffected.

This separately versioned, one-CPU saved-only analysis replaces that descriptor
with **inter-residue** mapped-heavy distances, excluding same-residue pairs and
the existing graph's1–2/1–3 neighbors. It does not attempt to complete the PQQ
bond graph or call these distances clashes. It retains the twelve closest pairs.

Also compare origin and each saved proposal against full-source heavy atoms
omitted from the exact physical source map. Restrict the queried donor atoms to
mapped standard Asp/Glu/Asn sidechain heavy atoms; exclude backbone and the full
same residue, so no inter-residue peptide/covalent bond is interpreted as contact.
Report every role's admitted source context, not a score-selected residue. These
standard sidechains have no inter-residue covalent bonds in the declared model.
Require the already-measured mapped-origin/full-source coordinate-frame residual
to be <=1e-10Å; otherwise mark this additional metric unavailable. Preserve the
actual donor displacement and atom names with each nearest pair. A donor that
did not move remains explicitly zero displacement, not evidence for packing
relief. Neither a short distance nor a change proves an omitted interaction is
energetically important. Spectator metals are excluded from this heavy-contact
descriptor, not deleted from the source or energy model.

Make a compact48-row TSV of existing confidence, coordination, selection, motion
and recorded boundary/residual fields. No new molecular calls, coordinate
optimization, folds, thresholds, calibration or production changes. All numerical
execution remains Slurm-only. Root owns submission and interpretation.
