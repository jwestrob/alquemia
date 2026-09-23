# Uniform continuation repairs the tested abstention; five cases remain unstable

The completed uniform32 test has a useful discrimination result and an explicit
qualification failure. Both prescribed passes retain25/25 canonical calls and
3/3 consumed crystal calls under the frozen cold-protocol bands, and change the
four-fold development set from3correct/1abstention to4correct. However, five
canonical sources fail the predeclared numerical checks. The new reference
correctly remains unavailable; this is not a promoted or fully qualified scorer.

## Same fixed geometry pools, all cells included

| Group | Cold old-band calls | Pass2 old-band calls | Numerically qualified pools |
|---|---:|---:|---:|
| Designated canonical references |25/25 correct|25/25 correct|20/25|
| Consumed crystals |3/3 correct|3/3 correct|3/3|
| Four specified noncanonical folds |3correct,1abstention|4/4 correct|4/4|
| Total |31correct,1abstention|32/32 correct|27/32|

Mathematical and operational decisions agree in this table. These raw old-band
transfers are distinct from a calibrated new-protocol result. All32 sources,
all384 cells and both passes remain in the denominators. This is consumed
development evidence; the folds and canonical sources are correlated, and
these functional-class labels are not direct affinity measurements.

Q4W6G0 Ca-conditioned sample2 changes from R=−405463.62238655 to
−405471.61668658 model kcal/mol, restoring its expected Ca-supported call.
Its La row now selects adaptive_La instead of adaptive_Ca after the vacuum
initialization discrepancy settles. This is the full uniform-pool result,
not replacement of only a favorable cell. The other three specified folds
retain their correct calls. The second changed endpoint choice is canonical
P38539 La: adaptive_La→adaptive_Ca; its pool remains Ca-supported but fails
numerical qualification.

Descriptive raw canonical class gap is7.32538885kcal/mol cold and7.33655917
after pass2: essentially unchanged, not a broad improvement in score spread.
No new bands were assigned from those raw extrema. The frozen rule requires
all25 designated calibration members to be qualified; only20 are. Both new
mathematical and operational reference records therefore have null bands and
`unavailable_calibration_member` status. Crystals and folds never enter fitting.

## Five specific numerical failures

379/384 scalar cells satisfy |pass2−pass1|≤0.1kcal/mol. The failures are:

| Canonical source | Candidate / metal / medium | Pass2−pass1, kcal/mol |
|---|---|---:|
| C5B120 |adaptive_La /La /ALPB|−0.631867|
| BBL57595.1 |origin /La /vacuum|−0.637313|
| P15279 |adaptive_Ca /La /vacuum|−0.353550|
| P38539 |adaptive_Ca /La /vacuum|+4.110974|
| Q4W6G0 |origin /La /vacuum|+0.121636|

Four of96 same-geometry Ca−La contrasts exceed0.2kcal/mol. Three of32 pooled
contrasts exceed0.2: C5B120(+0.645138), P15279(−0.355137),
P38539(+4.112404). BBL57595.1 and Q4W6G0 pass the pooled check but retain
their failed required component checks. Nothing was dropped because another
candidate or the final classification looked acceptable.

All768 logical outputs declare normal SCF convergence and satisfy the printed
energy-change tolerance. **None satisfies the printed MAX-density tolerance;
only16/768 satisfy RMS-density tolerance.** Actual residuals, electronic
components and output receipts remain available. The small eight-call success
therefore does not make two self-continuations a general density-convergence,
force-consistency or global-minimum guarantee. The frozen second pass is used
even when it raises energy; no third pass or favorable-stage selection occurred.

## Actual execution and cost

All384 exact cold GBW+xtbw pairs were available. Six cells reused exact completed
two-pass outputs:1H4I origin across both metals/media, Q4Ca2 La/adaptive_Ca/vacuum,
and P385La2 La/adaptive_La/vacuum. Reuse required the same coordinates, state,
medium, parameter export and both cold seed hashes. Those12 calls retain their
historical rank/receipt identity. Every other cell received both continuations.

Job1211168 completed **756 new scalar calls**, zero execution failures, with
actual XTBRESTART confirmed for every cell. Same nativeGFN2/ALPB or vacuum,
electronic300K, MaxIter500, native mixer, charge/protonation/water inventory and
fixed geometry. No new MACE, DFT, optimization, cold recomputation, new seed
strategy, source substitution or changed production input.

Actual allocation: **422s ×32CPUs =13,504 allocated core-seconds**,64GiB
requested, zeroGPU. Engine runtimes sum7,356.551s over756 one-rank calls,
overlapping across32workers. The difference includes seed copying, manifest/
receipt/parameter validation, collection and calibration inside the allocation;
it is part of measured cost. Peak batch RSS reported5,876,316KiB. Local
preparation/testing/report time and historical reused calls are separate.
This is an added numerical protocol cost, not a demonstrated production speedup.

Four final real-artifact tests pass, zero skips in11.088s: exact population/
self-seed reuse, finite manifests and geometry pins, archived canonical
calibration replay, and recomputed actual384-cell/32-pool/25-reference algebra.
Prelaunch explicitly left the new scientific result untested. Inventory v1's
origin-audit metadata lookup issue was repaired from the actual source outputs;
v1 remains preserved and no molecule was rerun to fix that lookup.

## Delivery and recommendation

Actual outputs under
`workspaces/precision_pool_continuation_20260923/run_v1/`:
`COMPARISON.json`, `REFERENCE.json`, `COSTS.json`, both stage manifests and
collections. The comparison preserves raw matrices, choices, every component
shift and each failed gate. The reference is explicitly unavailable.

Retain the unchanged existing scorer. Preserve this continuation route as a
promising numerical candidate: it produced a real abstention repair and removed
the diagnosed tiny-geometry discontinuities, but the remaining canonical
instability needs a distinct decision before integration. Do not launch more
identical passes or a production rescore from this report. This bounded test is
complete, and no further calls are scheduled by this branch.
