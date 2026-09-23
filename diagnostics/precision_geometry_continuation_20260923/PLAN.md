# Exact anomaly cells: two fixed native self-continuations

## Question and authorization

Do the existing two-pass native GFN2 self-continuations settle the electronic
solutions that change by approximately 10–11 kcal/mol under tiny geometry
differences? The completed exact geometry × rank test (cff57c9) excluded rank
count as the cause in these two cells. Prior four-pool continuation did not
include these exact geometries.

Jacob's standing overnight authorization applies. Root explicitly authorized:
“assess whether previously implemented fixed two-pass native continuation
settles THESE geometry-sensitive solutions ... four rank1 targets (two
sources×old/new geometry) ... 2passes each=8newGFN scalar calls ... Uniform
finalstage regardlessenergy ... do not cross-seed or widen beyond this
prescribed8-call diagnostic.” Following inventory, root confirmed: “proceed
exact8calls/fourpairedrank1seeds/twofixedpasses ... CPUonlyfourworkers”.

## Fixed population and method

Use both old/new geometries for each exact rank-one cold result in
`workspaces/precision_rank_geometry_20260923/run_v1/COLLECTION.json`:

- Q4W6G0, Ca-conditioned sample2, La at adaptive_Ca, vacuum.
- P38539, La-conditioned sample2, La at adaptive_La, vacuum.

Four cells; all charge −1, multiplicity1. Preserve exact XYZ bytes and existing
ORCA6.1.1/native GFN2 recipe, native mixer, electronic300K, MaxIter500,
parameter exports and convergence settings. The only initialization change is
the previously verified matched GBW+xtbw AutoStart route. Each stage1 starts
from its own exact rank-one cold seed pair; stage2 starts only from its own
confirmed stage1 pair. No cross-geometry seed, stage selection by energy,
TightSCF, third continuation, new optimization, MACE, DFT or gradient call.

Require actual `INITIAL GUESS: XTBRESTART`, native mixer, same parameter
export/electron count/charge/multiplicity, successful receipts and preserved
seed-before/seed-after records. Unavailable/failed stage1 cannot seed stage2;
retain all four cells and null missing values. There is no fallback.

## Frozen comparison

Always report stage2, preserving cold and stage1 values. The necessary settling
criterion is |E(stage2)−E(stage1)| ≤0.1 kcal/mol per cell, matching the existing
component policy. Report signed new−old geometry work for cold, stage1 and
stage2, and whether the stage2 geometry difference is within0.1 kcal/mol.
This latter check measures continuity at these near-identical coordinates;
it is not evidence of a global electronic minimum or a new affinity score.
Report both gates separately. No full pool, class assignment or refitting.

## Finite execution and cost

At most eight new scalar calls: four concurrent one-rank workers, two dependent
stages, four allocated CPUs/16GiB, CPU-only GPU-partition route on the existing
explicit host. Existing native runner and seed/receipt parser are reused; new
thin adapter restricts scope and one-rank execution. Earlier exact2×2 allocation
finished in19s, with individual calls in the seconds range; this should be a
short CPU diagnostic. Report actual allocation and call times, including
validation/staging, without conflating historical reused work.

No production/default/reference/historical-result changes. This plan is frozen
before the eight molecular calls. Stop after the fixed two stages and report
failure honestly if continuity or settling fails.
