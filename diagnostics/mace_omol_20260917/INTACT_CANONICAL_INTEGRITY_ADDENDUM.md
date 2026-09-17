# 1KB0 whole-chain preparation invalidated before its MACE evaluation

2026-09-17: a direct geometry audit found two inferred peptide bonds across
unresolved structure: C511–N513 is 4.750089 Angstrom and C573–N579 is
18.969263 Angstrom. The archived protonation manifest reported no missing
residues, and OpenMM's template matching accepted these inferred connections.
Template matching and preserved heavy coordinates were therefore insufficient
to validate this whole-chain preparation. The terminal OXT addition did not
repair these internal gaps.

This was discovered while job1200819 had started only its first calibration
protein, before any 1KB0 endpoint. All 25 calibration structures and both earlier
crystal transfers pass the same peptide-distance inspection. No 1KB0 score was
used to identify or choose this exclusion.

Keep the existing job, manifests and outputs intact. Its four 1KB0 tasks may
produce raw diagnostic outputs; they are **not valid benchmark scores**. The
report must enforce the independent preparation audit, set the 1KB0 whole-chain
score and class unavailable, and retain it in the three-transfer denominator.
The other 25 new cases and eight valid earlier endpoint reuses remain useful.
No label, charge, threshold, water state or source geometry is changed.

The audit checks every recorded inter-residue backbone C–N bond against its
actual heavy-atom coordinates, accepting 1.0–1.8 Angstrom. This broad integrity
range is declared after observing the obvious structural gaps, before the
affected prediction; it is not a fitted predictive criterion. Record all
distances, source identities and preparation hashes. Unsupported preparation
overrides numerical success in reporting. Do not use the original frozen
reporter without this new audit: use a new frozen reporting source.

The old fixed-core 1KB0 benchmark and baseline remain unchanged. This finding
concerns its new whole-chain representation, not evidence against the already
audited canonical core. A complete source/model or a separately declared and
validated reconstruction is needed for this whole-chain transfer. Do not
silently bridge the gaps, cap them, or substitute a different structure.
