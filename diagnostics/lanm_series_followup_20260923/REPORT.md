# LanM La/Dy: cheap relative-selectivity signal, with site dependence

**The completed model supports stronger relative La preference in Hans than
Mex at EF1 and EF3; EF2 is nearly neutral and slightly opposite.** This is a
useful within-series electronic signal on real structures. It is not a
quantitative affinity result or three independently labeled site predictions.

All **12MACE +48nativeGFN2** calls succeeded, including24 positively confirmed
self-continuations. All12 final composite endpoints are available. No DFT,
optimization, new waters, fitted weighting or production/default change occurred.

## Actual ordered result

Each endpoint is native OMOL vacuum plus native GFN2(ALPB−vacuum). The balanced
exchange is

`D_i = [E_Hans,i(Dy)−E_Hans,i(La)] − [E_Mex,i(Dy)−E_Mex,i(La)]`.

PositiveD means the model makes Hans more La-selective than Mex relative toDy.
Aqueous metal references and element-dependent atomic offsets cancel. Each
protein's atoms/water inventory appears unchanged on both sides of the cycle;
there is no missing aquo energy replaced by zero. Metal-dependent relaxation,
protonation and population effects remain omitted.

| Ordered site | Native MACE component | Solvent component | CompositeD |
|---|---:|---:|---:|
|EF1|+5.991489|−1.125121|**+4.866368**|
|EF2|−0.003447|−0.260992|**−0.264440**|
|EF3|+13.022013|−0.522897|**+12.499116**|

All values are modelkcal/mol. Full unrounded endpoint values, component energies,
native state audits and receipts are in the immutable final_collection.json.
The solvent term modestly weakens Hans's relativeLa preference; the main signal
comes from native MACE. We do not choose EF1/EF3 and discard EF2, sum the sites,
fit weights, interpret these magnitudes as measured free energies, or apply
La/Ca bands. There is **one consumed protein-level comparison**, with an ordered
three-site description and no independent site-affinity labels.

## What experiment supports—and does not support—the comparison

The [primary Hans/Mex study](https://www.nature.com/articles/s41586-023-05945-5)
reports Hans main apparentLa affinity68pM and Dy2.6nM atpH5; Dy gives partial
folding and heterogeneous site responses. Its Fig1e compares Hans with previously
published Mex pH5 CD apparent-affinity results from
[Deblonde etal.2020](https://doi.org/10.1021/acs.inorgchem.0c01303). The primary
conclusion is substantially stronger Hans light/heavy discrimination.

This is the same pH and apparent-affinity observable, **not one contemporaneous
assay batch**. The exact Mex La/Dy numerical values and uncertainties were not
recovered in this turn; supplementary retrieval failed. We use the primary
qualitative direction, not an approximated magnitude target. CD affinities
include folding/cooperativity and Hans dimerization; isolated electronic sites
do not reproduce that entire process. The supplied positiveD values are much
too incomplete to convert to Kd or a protein-level selectivity factor.

## Actual source/state and capability checks

Hans uses the preserved50atom pH5 peptide-amide-v3 cores from La-bound8DQ2 chainA,
sites201–203, no first-shell waters. Mex uses Nd-bound8FNS chainA, sites201–203,
44atoms/site and two deposited waters each:326/345,322/328,320/331. Source-heavy
positions are unchanged; La/Dy nuclear pairs and formal charge−1 are exact.
Their different crystal-conditioning states are an important unresolved source
of apparent selectivity, motivating structural transfer rather than promotion.

The first Mex preparation correctly failed on a remote Ser52 alternate-conformer
mismatch. A separate version applies the existing occupancy-aware selection
before protonation and preserves every deposited heavy coordinate exactly.
Both attempts remain archived; no geometry was adjusted to improve a score.

Actual OMOL input batches verifyLa multiplicity1 andDy multiplicity6, charge−1,
the unmodified checkpoint/head and expected atom count. This confirms the chosen
total-spin input, not localized f-electron physics or a verified ground state.
Native GFN2's explicit3-valence-electron5d/6s/6p model treats4f electrons in core:
both solver inputs use effective multiplicity1. Actual output counts are150
valence electrons for each Hans endpoint and136 for each Mex endpoint, with
charge closure and matching parameter exports. These physical and effective
states are recorded separately; production's Ca/La singlet guards were untouched.

## Numerical sensitivity

Each native cell used the same300K, native mixer, MaxIter500 settings for an
initialNoAutostart evaluation and one same-cell continuation. MatchingGBW+xtbw
seeds activated all24 actual `XTBRESTART` initial guesses. Continuations took
3–10SCF cycles. Reported energies always come from the continuation, even where
it is higher; no favorable solution was selected.

| Site | D from initial native cells | D from continued cells | Change |
|---|---:|---:|---:|
|EF1|+4.874252|+4.866368|−0.007885|
|EF2|−0.270593|−0.264440|+0.006154|
|EF3|+12.492300|+12.499116|+0.006816|

The largest absolute individual continuation change is0.008359kcal/mol. These
checks leave the qualitative ordered pattern unchanged. They do not establish
global SCF-branch convergence. All24 continuations still exceed the generic
printed maximum/RMS density tolerances, as seen in the earlier native-GFN
diagnostic; native normal convergence and empirical energy stability are recorded
separately. No extra continuation, threshold adjustment or failed-value fallback
was run.

## Cost, coverage and artifacts

| Job | Calls | Elapsed | CPUs | Allocated core-s |
|---|---|---:|---:|---:|
|1210308|1MACE capability|9s|32|288|
|1210309|1GFN initial capability|7s|8|56|
|1210311|1GFN continuation capability|7s|8|56|
|1210314|11MACE remaining|12s|32|384|
|1210315|23GFN initial remaining|44s|64|2816|
|1210319|23GFN continuation remaining|45s|64|2880|
|**Total**|**12MACE/48GFN**| | |**6480**|

MACE requested oneH200 for9+12s: **21requestedGPU-seconds**. Scheduler GPU TRES
fields were absent; both runtime receipts identify theH200. Native jobs requested
noGPU. Peak MACE allocatedVRAM898448384bytes and hostRSS1787856KiB; largest native
batchRSS1850056KiB. Local preparation/reporting is additional unmetered work.
There were zero molecular failures. One rejected source-preparation attempt is
preserved. The capability subset was reused within the totals.

Seven focused real-artifact tests pass in2.347s, zero skips. Tests and actual
scientific execution are distinct: the latter additionally confirms physicalDy
batch6, effective nativeDy singlet/electron count,24 native restarts, source
closure and complete paired results. No test fabricates successful molecular data.

- Result: `workspaces/lanm_series_followup_20260923/prepared_v2/final_collection.json`
- Costs: same directory,`COSTS.json`
- Compact vector/sensitivity: `RESULT.json`; exact pins:`ARTIFACTS.json`
- Method/authorization:`PLAN.md`,`EXECUTION_AGREEMENT.md`; operations:`COMMANDS.md`
- Last prior LanM test/history:`READINESS.md`

**Recommendation:** pursue the specifically authorized Dy-conditioned Hans
structural transfer before claiming transferable within-series discrimination.
It can test whether the signal survives a real change in conditioning geometry.
Keep the current result intact and keep the scanner default unchanged.
