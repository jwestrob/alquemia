# Saved Asn contact changes after bounded relaxation

Use the actual terminal remaining-stage collection, which retains the completed
4831 scout, and all ten prepared sources. Compare exactly four saved geometries
per source: original origin, fixed amide flip, Ca-optimized flipped proposal and
La-optimized flipped proposal. Both exchanged-metal cells must exist and have
identical physical coordinates for each measured geometry. Missing or failed
candidates stay unavailable, with original optimizer/boundary/failure evidence;
there is no substitute geometry or successful-source subset.

Track named Asn290 OD1/ND2–metal distances and the exact nonbonded pairs
HD22–Gly292 backbone H and HD21–Gly292 backbone H. Also track each amide H versus
Asn290's own backbone H, so the earlier peer contacts remain distinguishable.
Keep actual source identities/indices and before/after distances for each pair.
These named records, not a generic nearest-contact value, determine whether the
previously reported1.38Å target contact became longer or shorter.

Separately report the shortest included H-involving nonbonded pair, shortest H–H
pair, and shortest pair that shortened relative to the fixed flip, among the
four amide moving atoms versus the other included physical atoms. Reuse the
source graph and both amide-H parents. Exclude synthetic caps, the moving group
against itself and directly bonded pairs; do not silently exclude graph-distance-
two pairs or call every short/shortened distance a clash. This adds no admission
threshold. The existing final full-source gate and optimizer boundary are retained.

Reuse the pinned geometry helper from actual4626 for XYZ/key handling and require
origin/fixed-flip OD1/ND2 distances to reproduce its saved values. All coordinates,
state, source graphs and 204-atom identities are immutable; no geometry is built,
no model is imported, and no energy or optimization runs. Distances and comparisons
execute only on one Slurm CPU. Metadata freeze occurs after the actual collection
and execution receipt exist, preserving partial terminal outcomes if any.

Results include all ten source rows, per-candidate missing/boundary status, full
named pair evidence and a40-row CSV. A distance increase is descriptive contact
relief only; it is not an energy attribution, unconstrained minimum, equilibrium
rotamer or validation of a biological class.

## Real-scout interface check

Before the final collection exists, the same analyzer can bind actual4831's
completed Ca0 row and nine real unexecuted source rows. The one-CPU check requires
all four Ca0 geometries, paired coordinates and named source identities to pass;
all nine other rows must retain their two measured baselines and two unavailable
proposals, without synthetic coordinates or substituted results. Baseline donor
distances must reproduce actual4626. A successful check reports explicit partial
candidate coverage alongside passed interface checks; it does not claim that the
remaining searches ran. The final all-ten analysis remains a separate immutable
binding after the real remaining-stage terminal collection is available.
