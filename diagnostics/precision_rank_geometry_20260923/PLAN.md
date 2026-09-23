# Exact geometry by rank: two scalar sensitivity diagnoses

Root authorized exactly four missing calls on September23 under Jacob's
overnight permission, after reading the completed uniform-precision225 result.
This plan is frozen before new molecular outputs. It does not change a scorer,
reference, classification, geometry, chemical state, optimizer or starting guess.

Two consumed source cells were selected from existing results:

1. Q4W6G0 Ca-conditioned sample2, La at the Ca proposal, vacuum: the new
   correct→inconclusive case; observed pooled shift+7.98730437kcal/mol.
2. P38539 La-conditioned sample2, La at the La proposal, vacuum: the largest
   absolute pooled shift (−11.07964498kcal/mol), chosen by magnitude, not label.

Both cells have charge−1, multiplicity1. Reuse old-geometry/eight-rank and
new-geometry/one-rank outputs exactly. Complete each2x2 with old geometry at
one rank and new geometry at eight ranks: **four new nativeGFN2 scalar calls**.
Every input/XYZ is copied from its corresponding successful actual cell.
Same ORCA executable, nativeGFN2, electronic300K, MaxIter500, mixer,
convergence tolerances, maxcore2000, vacuum and NoAutostart. No GBW/xtbw seeds,
continuation, optimization, MACE, DFT, threshold changes or automatic retries.
Retain failed calls explicitly; never choose the lower or favorable result.

The read-only inventory inspected616 existing tasks across the rank28 panel,
rank pilot, two native continuation stages and three native-force stages.
None had an exact matching geometry/state for the two targets. Thus existing
continuation does not answer these particular cells. Exact actual source/recipe,
output/receipt/parameter pins are in
`workspaces/precision_rank_geometry_20260923/INVENTORY_v1.json`.

Before seeing results, declare the same0.1kcal/mol component gate for an
exact-geometry rank comparison. Report E(new,r)−E(old,r) separately at each
rank, and E(g,1)−E(g,8) at each geometry; also report their interaction. Values
remain visible on failed gates. These four energies can diagnose sensitivity,
not establish a unique electronic ground state or an improved classifier.
Retain actual SCF cycles, electronic/state audits, parameter bytes, timing and
all normal/failure receipts. No updated pool or calibration is produced.

One CPU-only allocation on the existing H200 host:18 CPUs/64GiB, no GPU.
The existing runner requires a uniform rank count per manifest, so prepare
two disjoint contained manifests, two tasks each: two one-rank workers plus
two eight-rank workers run concurrently (18 ranks total). Each has separate
locks and events. Both preserve actual18CPU allocation metadata; executor
accounting overlaps and must not be summed. Authoritative allocation cost is
one SLURM allocation wall time×18, with no duplicated cost or molecular calls.
Existing individual target runtimes are19.095/29.052/61.580/29.703s; this is a
small minutes-scale diagnostic, not a new full panel.
