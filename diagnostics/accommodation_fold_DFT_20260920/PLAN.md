# Preserved DFT baseline on the known-reference fold challenge

Frozen 2026-09-20 before new DFT results. Parent assignment under Jacob's
standing overnight accuracy authorization explicitly requests this experiment.

## Question and inputs

Does the preserved DFT discriminator tolerate the same structural conditioning
challenge as the MACE variants? Use **all 250 frozen sources** for the existing
25 experimentally associated PQQ reference proteins. These are consumed
calibration proteins and correlated structural repeats, not new blind labels.

Input: `workspaces/accommodation_goal_20260920/folds_v1/preparation_reconciled_v1.json`.
Of 233 supported preparations, reuse the 25 canonical original DFT pairs after
recipe/state/coordinate/receipt checks. Execute **416 new endpoints** for all
208 supported noncanonical sources (98 La-conditioned,110 Ca-conditioned).
Keep all 17 unsupported sources explicit in the 250/225 denominators. No
repair/replacement, new protonation, optimization, context DFT or changed water.

Copy the exact already-emitted original-core Ca/La inputs and XYZ files without
editing their scientific contents. Retain the existing native ORCA6.1.1
`r2SCAN-3c NoAutostart CPCM(Water) DefGrid3` recipe, default NormalSCF,
`maxcore8000`, native basis/ECP, original charge/multiplicity and dry PQQ core.
Only allocation-specific PAL16 is added through the existing pinned renderer.

## Frozen decisions and aggregation

Use the authoritative released record
`diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json`:

- `R=(E_Ca-E_La)*627.509474`, kcal/mol; larger is more La-like.
- Ca-supported when `R <= -405404.9238187139`.
- La-supported when `R >= -405396.32088671124`.
- Between these bands: inconclusive. No threshold refitting.
- Optional aquo-gauge S is display only, using the existing pinned aquo gap;
  it is not a new absolute-affinity reference or decision threshold.

Report all single sources, canonical25 replay, primary225, La100 and Ca125.
For each protein: median of **all four noncanonical La sources**; median of
**all five Ca sources**; and the equal mean of those two arm medians. Any
missing member makes that arm, and any dependent balanced summary,
unavailable. Do not select successful members. These are robustness
descriptors under frozen single-source bands, not thermal free energies or
a newly calibrated ensemble classifier.

Later compare side by side with the matching MACE comparison on the exact
same preparation manifest. Preserve core/context and native/composite MACE
results separately; no flexible fitting or claiming improvements from25/25.

## Execution and real prior cost

Four disjoint manifests of104 endpoints; four64-CPU standard,memory jobs,
each with four16-rank ORCA endpoints concurrently. Preserve maxcore8000 and
request600GiB per node to cover declared per-rank memory plus overhead.
Use the existing64-CPU node exclusions and MPI/thread policy. No custom
time/compute cap; scheduler ceilings still apply. No new executor or default.

All50 original compatible endpoint receipts are available: total wall
17,822.282511 s; median349.3619455 s; range153.731546–555.134219 s;
267,334.237665 endpoint-rank-seconds at15 ranks on node-344-8t-1. These are
matched scientific-core/method costs, **not matched-hardware predictions**
for the new16-rank/64-CPU batches. Preserve individual measurements in the
prepared manifest; measure new throughput rather than citing queue estimates.

Keep every raw output, unrounded endpoint energy and execution receipt.
Collectors run even if an executor exits nonzero, record partial/failures,
and never substitute cached baseline scores for failed new endpoints.
No adaptive scientific rescues are authorized in this finite manifest.
