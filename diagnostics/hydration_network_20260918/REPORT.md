# Water preparation improves the original La/Ca discriminator

**Correcting water hydrogens repairs the alpha-lactalbumin/GGR ordering in both
alpha structures, using the original scoring cores and DFT recipe.** MACE now has
a demonstrated role as a cheap geometry proposer for this consumed development
case. Production/default and PQQ inputs/results remain unchanged.

## Main result: isolate preparation from the energy model

The original40/43-atom alpha cores retain exactly their protein, caps, metal,
water oxygens, atom order, water count and charge. Only water H coordinates are
transferred from MACE optimization in source-defined expanded hydrogen-bond
networks. The scoring calculation uses the original native r2SCAN-3c/CPCM/DefGrid3
SP recipe and SCF policy; La0/Ca−1 singlets. It does not add environmental energy,
change a threshold, choose a favorable deletion or score extra protein atoms.

`R=E_Ca−E_La`. The table is `(R_alpha−R_GGR)*627.509474`, kcal/mol; positive is
the expected direction. No aquo offset or absolute classification threshold is
needed. All three GGR cores are water-free and their original compatible results
are reused after verifying inputs, outputs and receipts.

| Alpha structure | GGR structure | Original | Prepared water H |
|---|---|---:|---:|
| 1F6S | 1GLG | −14.559 | +10.924 |
| 6IP9 | 1GLG | −17.543 | +14.573 |
| 1F6S | 2FW0 | −24.828 | +0.655 |
| 6IP9 | 2FW0 | −27.812 | +4.304 |
| 1F6S | 2FVY | −23.648 | +1.835 |
| 6IP9 | 2FVY | −26.631 | +5.484 |

Thus0/6→6/6 structural comparisons, representing **one consumed biological
comparison between two protein groups**, not six independent successes. The
primary predeclared comparison was to1GLG; the other two previously inspected
GGR structures are a subsequent read-only robustness comparison. The weakest
margin is only0.655kcal/mol. This is a useful development improvement, not broad
validation, an experimental affinity magnitude or physiological occupancy.
Alpha's condition/construct qualification and the difference in assays remain;
see the [prior evidence discussion](../ggr_mechanism_plan_20260915/REPORT.md).

Unrounded results and receipts: [RESULT.json](RESULT.json),
[GGR structural comparison](GGR_REPLICATE_RESULT.json).
New scoring protocol: `amide_v3_native_r2scan3c_context_prepared_water_H_v1`.

## What was wrong, and what MACE contributes

The preceding pilot already fixed stretched water O−H bonds; existing amide-v3
already fixed the coordinating-backbone defect. Those were retained. Here,
source-neighbor selection restores omitted peptide/sidechain hydrogen-bond
partners and one/two outer waters in70/76-atom preparation contexts. Protein
fragment composition is common across the two alpha structural replicates.
Only the inner-water orientations move; all oxygens/protein/caps stay fixed.
Context endpoint charges are La−1/Ca−2; original-core scoring retains La0/Ca−1.
This is contextual geometry preparation, not a subtractive environmental energy.

Generated H atoms in the problematic A211/A310 waters initially approach the
metal to1.852/1.880A. A geometry-defined outward start lowers both DFT energies,
especiallyLa. Restoring neighbors plus tighterSCF shifts initialR only1.05/1.53
kcal/mol, whereas orientation shifts it24.82/31.16. This identifies a substantial
preparation mechanism; it does not uniquely decompose the electrostatic cause.

Native MACE-OMOL reproduced8/8 local DFT energy-change signs and4/4 starting-seed
rankings; median water rotational-gradient cosine0.991329 and local-change
MAE0.424043kcal/mol. The eight exact rigid-water MACE searches all converged,
with source/radial minima differing by at most1.1e-5kcal/mol per endpoint.
No retraining, charge masking, model selection by label or numerical DFT gradient.

Four native DFT energy/gradient checks confirm every selected MACE configuration
has lower DFT energy than either original start. DFT rotational gradients remain
3.47–5.42kcal/mol/radian: these are improved preparations, **not DFT stationary
points**. The original long native-DFT searches1201824/1201825 remain running as
comparators. Their small frozen-coordinate drift and absent final.engrad are
explicitly handled; no convergence/geometry criterion is silently relaxed.
Equivalence to fully optimized DFT minima is not established or required to state
the original-core preparation result above.

## Actual cost and checks

Practical preparation/scoring route for both structures:

- MACE:96s, oneGPU,16CPUs; eight searches/136 energy-force objective evaluations.
- Original-core DFT:118s on64CPUs, four16-rank endpoints, no retries.
- Combined allocation:9,088CPU-core-seconds and96GPU-seconds. Preparation/tests
  ran separately and were not individually metered. Per-endpoint DFT wall times
  are106.06–115.47s, versus74.12–99.16s in the archived matching alpha endpoints;
  this is an indicative comparison, not a controlled hardware speed benchmark.

One-time development costs are separate: proposal validation137GPU-s/2192core-s;
expanded-core DFT adjudication402s/25728core-s; shared bulk-water reference25s/
1600core-s including a one-second MPI-startup failure. Ongoing native-DFT searches
are additional development cost, not included in the proposed routine route.
[Completed scheduler receipts](COMPLETED_ACCOUNTING.txt) retain exact accounting.

**40 real-fixture tests pass, zero skips:**13 hydration/projection/parser/transfer,
four prior water regressions,23 archived-baseline and preparation tests. Released
PQQ energies and bands pass regression without new quantum calculations.
All eight new protein DFT single points/gradients completed normally. The first
expanded-core collection encountered a parser ellipsis in gradient progress
text; the numeric regex was fixed and recollection succeeded, with no QM rerun.
Old failed collection and execution receipts remain. A missing collector import
was restored from its pinned existing dependency, with no quantum input change.

## What remains

Water identity is now modeled through source atoms and local hydrogen-bond
partners; orientation is metal-dependent under the same rules. Occupancy is
still a separate problem. The bulk-water reference is implemented, but bound
water entropy/ZPE shifts and non-electrostatic terms remain unavailable. No
occupancy probability or binding free energy is invented.

Next: collect native orientation comparators, then use the useful cheap proposal
route to evaluate the declared joint occupancy patterns. The36-task DFT-only
occupancy manifest remains unsubmitted and is superseded as an execution choice
pending that assessment. Extend the preparation rule to additional water-bearing
benchmark cases before considering default promotion. Keep the successful
baseline available and preserve dry PQQ fidelity.

Operations: [COMMANDS.md](COMMANDS.md). Implementation/result email was sent to
Jacob through the previously authorized relay; local acceptance is recorded.
Vault note: `2026-09-18_laca-water-networks-and-mace-proposals.md`.
