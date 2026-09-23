# Fixed context plus adaptive accommodation improves individual fold calls

The completed frozen225-source transfer corrects both wrong released PQQ calls.
It returns **205 correct, zero wrong, one inconclusive and19 unavailable**.
Against the recovered original-context adaptive method, it corrects one further
wrong call. This is a useful, localized improvement on consumed reference
structures; it does not establish new biological affinity validation.

The25 designated canonical references and three consumed crystal controls
remain correct under the separately frozen union/adaptive reference. No transfer
source entered calibration. Production, explicit DFT access, historical scores
and the unknown PLM cohort remain unchanged.

## Individual sources, with all225 retained

| Method, own frozen reference | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Released static context composite |203|2|2|18|
| Recovered original-context minimal adaptive |204|1|1|19|
| **Fixed-union minimal adaptive** |**205**|**0**|**1**|**19**|
| Fixed-union static |199|1|8|17|
| Standalone-xTB static composite |201|2|5|17|
| Preserved native DFT |199|3|6|17|

On205 sources available in both new and released methods,201correct/2wrong/
2inconclusive becomes204correct/0wrong/1inconclusive. The two errors corrected are
A0A3F2YLY8 Ca-sample1 and A0ACD6B9F2 Ca-sample4. The two corrected inconclusives
are A0A3 Ca-sample3 and A0AC La-sample4. **C5AXV8 La-sample3 becomes inconclusive**
(R=−405457.87754655955modelkcal/mol); it was already inconclusive under the
original-context adaptive method. On205 shared sources versus that method,
203correct/1wrong/1inconclusive becomes204/0/1.

A8R3S4 Ca-sample3 is available/correct here but unavailable in the released and
recovered original-adaptive ledgers. Its prior static-union recovery already
exists; this technical coverage gain is not evidence from new accommodation.
Two other previously scoreable structures are lost to optimizer stopping.

## Fold aggregation and spread

La4 and Ca5 mean the four noncanonical La-conditioned and five Ca-conditioned
input folds, irrespective of the protein's experimental class. Every required
member must be available. Balanced is the equal mean of complete arm medians;
triples are every three-of-four La subset, with no favorable subset selection.

| Strict descriptor | Released correct / missing | Union-adaptive correct / missing |
|---|---:|---:|
| La4,25 proteins |23 /2|22 /3|
| Ca5,25 proteins |21 /4|21 /4|
| Balanced,25 proteins |20 /5|19 /6|
| All100 correlated La triples |94 /6|91 /9|

Every available new aggregate is correct. Released aggregate decisions were
already correct; this experiment adds no aggregate accuracy and reduces coverage.
Across42 common complete conditioning arms,20 score ranges narrow and22 widen
versus released; median range9.98717→10.63428modelkcal/mol. Against original
adaptive,23 narrow/18 widen/1equal, median11.14084→10.63428. Better individual
calls are not a uniform reduction in structural score variation. Folds and
triples remain correlated samples of25 proteins, not independent observations.

## What ran and what failed

Four disjoint51-source shards supplied408 unchanged SLSQP searches.406 were
admitted; two reach the200-iteration limit despite valid physical geometry:

- Q92WY9 Ca-sample2 / La:2061 evaluations,620.961500s.
- Q60AR6 La-sample0 / La:2052 evaluations,484.805682s.

Their final geometries, forces and traces are preserved; no rescue enters this
version. Together with17 existing preparation exclusions they explain all19
unavailable sources.175 of408 endpoints reach a physical boundary. Median search
wall time is4.342458s and median optimizer evaluation count17, showing how strongly
these numerical tails affect practical cost. The separate precision pilot is
not imported into this frozen experiment.

Actual new calls:12,750 nativeMACE search evaluations,404 cross-metal MACE single
points and1,616 native GFN2 calculations. Every started MACE/GFN2 call succeeded;
optimizer failure is a separate status. All206 admitted source pools completed,
including four exact noncanonical pilot reuses.408 q0 force arrays and compatible
origin solvent components were reused. Zero new DFT, folds, protonation variants
or scientific retries. Both operational and mathematical pool scores are exactly
identical on this transfer; neither hides an alternative decision.

## Energy accounting and limitations

Both metals score the same finite geometry pool. E_M=E_OMOL,vac,M+
E_GFN2,ALPBwater,M−E_GFN2,vac,M; R=E_Ca−E_La. No mechanical/entropy scalar is added.
The canonical-only bands remain Ca_max=−405463.7090572145 and
La_min=−405456.46881754074, gap7.240239673759788modelkcal/mol. Larger R is more
La-like on this protocol's own scale; no absolute aquo reference or affinity is
claimed. Raw component work and all rows are retained in COMPARISON_v1.json.

Union membership was fixed without labels from ten saved folds per protein.
This exact preparation is not yet qualified for a future three-fold scanner.
Context composition/cavity, proposal forces and finite-pool selection all
contribute; these results do not isolate one hydrogen bond or environmental
mechanism. Native solvent derivatives remain unqualified from the separate
force experiment; this branch uses nativeMACE forces only.

## Cost, verification and recommendation

New transfer allocations total **513,768CPU-seconds and4,123 requestedGPU-seconds**.
All molecular jobs ran on node-224-2t-8gpu-1; the one-core final comparison used
node-128-512g-8gpu-1. Four concurrent local preparations additionally used
252.46/246.68/249.23/252.54s wall time and878.71 userCPU-seconds: see exact process
receipts in COSTS_v1.json rather than treating their sum as elapsed wall time.
Scheduler totals include in-job preparation, failed searches and collection;
nested eight-rank receipts are not added again. Full branch development including
prior pilot/canonical allocations totals583,099CPU-seconds/4,943GPU-seconds,
excluding earlier reused calculations and local setup/analysis. This is not a
matched fresh-production throughput benchmark or folding-inclusive cost.

All8 final real-fixture tests pass in15.161s, zero skips. They replay all225/75/100 results, independently recompute
both energy-selection policies and strict medians, retain missing members and
check exact corrected/abstained/optimizer-failure IDs. See the final test log and
[commands](TRANSFER225_COMMANDS.md). No fabricated final scientific output is used.

**Recommendation: pursue this challenger.** It corrects the remaining individual
fold error without sacrificing canonical PQQ fidelity. Before scanner promotion,
qualify a uniform numerical stopping revision and a deployable fixed-membership
preparation policy; keep failures visible and the released default accessible.
