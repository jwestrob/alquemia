# Uniform-precision transfer: complete, with a solvent-sensitive regression

**206 correct, zero wrong, two inconclusive,17 unavailable across the original225
folds.** All208 physically prepared sources scored. Compared with the original
union/adaptive result, two optimizer failures recover as correct calls and one
previously correct call becomes inconclusive. This improves coverage; it does
not establish numerical equivalence. Production and all historical ledgers are
unchanged. No result-triggered retry, threshold edit or additional source ran.

## Classification and coverage

| Frozen method | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Archived DFT |199|3|6|17|
| Released local-context composite |203|2|2|18|
| Original union/adaptive |205|0|1|19|
| Uniform-precision union/adaptive |206|0|2|17|

On the206 sources available under both union/adaptive policies:204 remain
correct, one remains inconclusive and one changes correct→inconclusive. The
two recovered cases are Q92WY9 Ca-conditioned sample2 and Q60AR6 La-conditioned
sample0, whose successful precision34 pools were reused exactly. C5AXV8
La-conditioned sample3 remains inconclusive. **Q4W6G0 Ca-conditioned sample2 is
the new inconclusive case**; the component audit below explains the observed
change without assigning an unproven unique cause.

The new candidate uses the already frozen canonical25 reference
`1f8470bbdf7a056098ca261aa10b72bdc940cf71ab860fbbc54090bda1e73250`:
Ca_max−405463.71230032056 and La_min−405456.3869114709 model kcal/mol.
Its prior25/25 canonical and3/3 crystal results remain preserved. No new
calibration occurred here. Every225 decision is also unchanged when this new
result is transferred through the old union/adaptive bands; the regression is
not created by recalibration. Mathematical and operational selections yield the
same225 decisions. Raw electronic contrasts are not binding free energies.

| Strict aggregate | Original union: correct / unavailable | Precision: correct / unavailable |
|---|---:|---:|
| La4 groups,25 denominator |22 /3|23 /2|
| Ca5 groups,25 denominator |21 /4|22 /3|
| Equal mean of complete arm medians,25 |19 /6|21 /4|
| Three-of-four La subsets,100 |91 /9|94 /6|

These aggregates have no wrong or inconclusive calls. They require all declared
members; missing structures remain missing. The100 triples and225 folds are
correlated structural samples, not independent biological observations. All
sources are consumed development. This supports improved completion and
structural robustness on this panel, not broad unseen-protein validation.
All17 original preparation exclusions remain. The independently repaired-H
MMOL1770 supplement is intentionally outside this ledger.

## Numerical and physical checks

202 fresh sources plus six exact reuses preserve the source geometry, chemistry,
union membership, q0 values/forces and paired four-angle selector. The single
optimizer change is the previously declared ftol1e−8 Hartree-equivalent.
The qualified one-rank scalar profile applies only to fresh solvent cells;
cached eight-rank receipts remain explicitly recorded.

- 404/404 new searches succeeded; maximum27 iterations and31 function evaluations.
- 414/414 comparable successful-old native proposal energies pass0.001kcal/mol.
  Maximum absolute shift0.00005742563kcal/mol; maximum mapped-atom movement
  between numerical policies0.00299849Å.
- 1624/1648 matched GFN2 medium cells pass0.1kcal/mol;198/206 pooled contrasts
  pass0.2kcal/mol. Combined source-level agreement is182/206, not206/206.
- Maximum absolute GFN2 component change44.51472448kcal/mol; maximum pooled
  shift11.07964498kcal/mol (P38539 La-conditioned sample2). The full JSON retains
  every failed gate, component, source and selection.
- 176/416 final endpoints have declared boundary flags (171/404 fresh).
  These are bounded candidates, not unconstrained minima or equilibrium populations.

All molecular attempts terminated successfully. Successful SCF convergence does
not eliminate the observed energy sensitivity. The rank28 check passed on its
actual controls; geometry and rank both differ in the fresh transfer comparison,
so that control alone cannot assign a unique cause to every transfer discrepancy.

## Q4W6G0: actual component and selection audit

Source `q4w6g0-pqq-la_model__conditioned_Ca__seed-1_sample-2`:

| Quantity | Original | Precision |
|---|---:|---:|
| PooledR, model kcal/mol |−405471.6096909145|−405463.62238654663|
| Call under either reference |Ca-supported|inconclusive|
| Ca-selected candidate |adaptive_La|adaptive_La|
| La-selected candidate |adaptive_La|adaptive_Ca|
| La-at-adaptiveCa vacuum,Eh |−258.884720999673|−258.867751100617|
| Vacuum SCF cycles /ranks |43 /8|29 /1|
| La-at-adaptiveCa ALPB,Eh |−259.115870722704|−259.115870844770|
| ALPB SCF cycles |29|29|

The La vacuum contribution changes+10.64877243kcal/mol while its matched ALPB
energy changes−0.00007660 and cross-MACE changes−0.00003953. Thus the solvent
transfer term becomes10.64884903kcal/mol more favorable, changes the La row
minimum, and shifts the pooled contrast+7.98730437kcal/mol. The relevant Ca
proposal coordinates differ by only2.82859e−6Å; the La proposal differs by
6.25995e−6Å. Their own native MACE energy differences are about1e−9kcal/mol.

Both vacuum outputs have normal termination, SCF convergence, return code0 and
452 effective electrons. All eight old/new input-recipe pairs and exported
GFN2 parameter pairs are byte-identical. The La Mulliken charge changes
0.87038089→0.77725458 in vacuum, compared with0.84291556→0.84291530 inALPB.
This establishes a difference in the converged electronic solutions, not a
missing atom, changed parameter set, parsing failure or threshold adjustment.
It does not distinguish geometry-triggered versus rank-triggered SCF sensitivity
or establish which electronic solution is appropriate. No extra calculations
were used to investigate this result. Exact outputs, receipts, charges,
SCF markers and all16 scalar cells are pinned in `Q4_COMPONENT_AUDIT_v1.json`.

## Actual execution cost

Four independent GPU shards used1H200/32CPU/200000MiB each; their dependent
native batches used32CPU/64GiB with32 one-rank workers. Jobs1211044–1211051 all
completed. Actual fresh calls:5935 search MACE,404 cross-MACE and1616 GFN2;
zero q0 orDFT calls and zero failed molecular attempts. The6339 optimizer
function evaluations include404 cached q0 evaluations; they must not be counted
again as fresh model calls. GPU summaries independently confirm the call counts.

**165472 allocated core-seconds and2312 requested GPU-seconds.** GPU allocation
durations were597/536/634/545s; CPU durations706/742/677/734s. Peak recorded
CUDA memory was5.682GB. The node ran disjoint allocations concurrently;
parallel elapsed latency is distinct from summed allocation. These totals
include in-job mapping validation, loading, collection and failure handling.
Local preparation/tests/reporting time was not separately metered; historical
six-pool reuse and calibration costs are excluded and preserved in their records.

For the same202 sources, optimizer function evaluations fall8962→6339 and
summed search wall time2339.414→1579.604s. Concurrent node load differs, so this
is actual cost evidence rather than a controlled hardware speed benchmark.

Nine final tests pass in9.965s, zero skips: real full225 source membership,
raw-matrix/selection replay, unchanged prior methods, frozen reference, qualified
rank profile, all75 strict groups and100 triples, exact six reuses, and physical
input containment. Pure replay tests are separate from the actual scientific
calculations above.

## Delivery and next decision

Full comparison: `workspaces/slsqp_precision_transfer_20260923/COMPARISON_v1.json`.
Cost receipts: `COSTS.json`; native response audit: `SEARCH_AUDIT_v1.json`;
Q4 output audit: `Q4_COMPONENT_AUDIT_v1.json` in the same workspace.
See [commands](COMMANDS.md) and the compact `ARTIFACTS.json` pins.

Retain this as a separately calibrated candidate with improved coverage and
the numerical limitations above. The Q4 regression remains visible. The next
useful diagnostic is to separate exact-geometry rank sensitivity from
geometry sensitivity; it must not pick favorable seeds or repair a test output
silently. Root separately coordinates fresh source integration and any further
scientific execution. No default or historical reference changed here.
