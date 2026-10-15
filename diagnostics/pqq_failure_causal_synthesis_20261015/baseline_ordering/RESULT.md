# One source pair explains the current accommodated class overlap

The completed frozen-charge common3 benchmark has one overlapping pair of
distinct sources: Ca-class ExaA Q9Z4J7 La-conditioned sample3 and La-class XoxF1
A0ACD6B9F2 Ca-conditioned sample4. Their ordering is reversed by
1.3704306245199405 kcal/mol. Both remain inconclusive under the fixed bands;
this overlap is not a confident wrong call. Operational and mathematical
selection policies give the same result.

The other current inconclusive sources create no additional class overlap.
In particular, the current Ca-conditioned125 cohort has an observed positive
class gap of1.444062333961483 kcal/mol despite four inconclusive memberships.
Thus inconclusiveness alone is not evidence of a malformed fold or an inverted
score ordering. It also reflects the deliberately unassigned reference interval.

| Existing cohort | Correct / wrong / inconclusive | Observed La-min minus Ca-max, kcal/mol |
|---|---:|---:|
| Expanded135 |129 /0 /6|−1.3704306245199405|
| Current Ca125 |121 /0 /4|+1.444062333961483|
| Two-Ca15 subset |13 /0 /2|+7.06399001320824|
| Combined260 memberships |250 /0 /10|−1.3704306245199405|

The combined260 memberships represent257 contexts. The A0AC Ca4 source appears
in both main cohorts, so the ten inconclusive memberships include the previously
reported nine distinct sources. These structures are consumed development data
and correlated replicas, not260 independent biological observations.

The unchanged canonical common3 reference gap is8.450323470286094 kcal/mol.
Observed cohort extrema are diagnostic only: no boundary was fitted, no decision
changed, and no narrower uncertainty interval is recommended from this analysis.
Origin-only results remain separate and have a larger negative class gap
(−6.857392679958139 kcal/mol) and two wrong memberships; the accommodated result
must not be attributed to origins.

## Consequence for the ongoing experiment

The already completed A0AC amide test addresses the low-scoring side of this
single overlap. The uniform transfer still has to establish what happens when
the same search is applied to every supported source. This saved diagnostic
does not establish the final common6 result, nor identify the physical cause
of the remaining near-band scores.

## Execution and reproduction

Saved-result job5317 completed on one Slurm CPU, with zero molecular calls or
new optimizations. Scheduler elapsed time rounds to0 seconds; measured script
wall time is0.1780953969573602 seconds. This is not a model-inference benchmark.

Result: `workspaces/pqq_failure_causal_synthesis_20261015/baseline_ordering_v1/execution/RESULT.json`.
Config SHA: `2d75f60d273eeb422d091ea4ae44769b7b7de204deb947444775d990fbb4f553`.
The immutable config pins the actual frozen-charge collection, current ledger,
script and declared cohorts. Exact submission and completion wake commands are
in the neighboring `SUBMISSION.json`. Outputs are immutable; a rerun needs a new
output directory and configuration, and no repeat is currently needed.
