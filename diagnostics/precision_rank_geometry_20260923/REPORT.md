# Exact2x2: the large energy changes follow geometry, not rank count

All four new native-GFN2 calls succeeded. At either fixed geometry, changing
one→eight MPI ranks changes the energy by at most6.43e−10kcal/mol. Both large
vacuum energy changes reproduce at both rank counts, including identical SCF
cycle counts. This removes rank count as the explanation for these two tested
discrepancies. It supports retaining the qualified one-rank execution profile;
it does not qualify an electronic ground state or change any classifier score.

| Source / La vacuum cell | Geometry | One-rank energy,Eh | Eight-rank energy,Eh | SCF cycles, both |
|---|---|---:|---:|---:|
| Q4W6G0 Ca2 /adaptive_Ca |old|−258.884720999672|−258.884720999673|43|
| Q4W6G0 Ca2 /adaptive_Ca |new|−258.867751100617|−258.867751100616|29|
| P38539 La2 /adaptive_La |old|−274.463644536509|−274.463644536510|39|
| P38539 La2 /adaptive_La |new|−274.481301114709|−274.481301114709|61|

Q4 new−old geometry work is+10.64877243kcal/mol at both ranks; P38539 work is
−11.07967010kcal/mol at both ranks. Maximum rank/geometry interaction is
1.25e−9kcal/mol. All four exact-geometry rank comparisons pass the frozen
0.1kcal/mol component gate.

Q4's Ca proposal moves only2.83e−6Å between old and new numerical policies;
the complete comparison already established negligible native-MACE and ALPB
changes. The vacuum outputs have different reported electronic charges while
both claim SCF convergence. The controlled2x2 therefore identifies sensitivity
to the near-identical geometry inputs, consistent with finding different SCF
solutions. It does not establish which solution is appropriate, or a physical
energy discontinuity. No energy was selected because it was lower or favorable.

## Exact scope and provenance

The fixed sources were Q4W6G0 Ca-conditioned sample2 (new abstention) and
P38539 La-conditioned sample2 (largest absolute pooled-contrast change).
Both cells are La/vacuum, charge−1, multiplicity1. Each new input and XYZ was
copied byte-for-byte from its actual source calculation. Same ORCA executable,
nativeGFN2,300K,MaxIter500,mixer,convergence,NoAutostart and parameter exports.
Only rank count changes within each geometry pair. Four old8/new1 cells were
reused; four missing old1/new8 cells were computed. Zero MACE, optimization,
DFT, new seed, retry, threshold change or updated pool.

The prior616-task inventory found no matching exact geometry/state in the
rank28 panel, rank pilot, previous pool continuation or native-force studies.
Thus this diagnosis required the four new calls; earlier continuation did not
already answer these two anomalies. Every source, output, execution receipt,
parameter audit and actual SCF cycle count is pinned in the collection.

## Cost and validation

CPU-only job1211105 completed in19s on18 allocated CPUs/64GiB, zero GPUs:
**342 allocated core-seconds**. Two contained uniform-rank runner manifests
executed concurrently: two one-rank plus two eight-rank workers. Their executor
wall accounting overlaps; it is intentionally not summed. The SLURM allocation
is counted once. This total includes setup, audit and collection; local
preparation/testing/report time and four historical reused cells are separate.

Four real-artifact tests pass, zero skips in the final run: exact input copies,
finite18-rank dispatch and containment, prior reuse inventory, and completed
matrix/sign/interaction algebra. Before execution the actual-matrix test was
explicitly skipped; no synthetic successful scientific output was used.

Actual collection:
`workspaces/precision_rank_geometry_20260923/run_v1/COLLECTION.json`.
This is a numerical diagnosis, not a recalculated or rescued classifier result.
The full225 historical comparison, its abstention and failed numerical gates
remain intact. Root separately coordinates any continuation test on these
specific electronic solutions.
