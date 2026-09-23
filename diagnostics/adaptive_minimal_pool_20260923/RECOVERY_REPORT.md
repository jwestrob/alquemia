# Completed minimal adaptive comparison, including recovered coverage

The three-geometry adaptive candidate retains its individual-fold improvement
while recovering two sources previously excluded by an obsolete proposal
dependency. On the same206 available folds, released scoring gives202 correct,
2 wrong and2 inconclusive; the candidate gives204 correct,1 wrong and1 inconclusive.
These are consumed structural repeats of25 reference proteins, not206 independent
biological observations. Aggregate reference decisions were already correct.

| All225 sources | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
|Released static|203|2|2|18|
|Minimal adaptive, archived replay|202|1|1|21|
|Minimal adaptive, completed origin recovery|204|1|1|19|

P12293 Ca-sample4 and Q60AR6 Ca-sample1 now both score correctly Ca-like. Their
static origins were already correct: these two calculations repair coverage,
not two classification errors. The original missing records remain beside the
new completed results. The separate Q60AR6 La-sample0 search failure remains
unavailable, as do17 preparation failures and one missing native origin.

The unchanged canonical25 reference also retains all three previously consumed
crystal checks. Strict aggregate completeness is preserved: La4 is22 correct/3
unavailable, Ca5 is21/4, balanced is19/6, and the100 prescribed triples are91/9.
No threshold was fitted to the recovered folds. Removing the old terminal-only
candidates changes11 raw scores but no available decision in the archived225.

## What ran and what this establishes

The [two-source molecular report](../adaptive_origin_recovery_20260923/REPORT.md)
records67 actual MACE energy/force calls and16 native GFN2 calls, costing5152
allocated core-seconds and49 GPU-seconds. All completed. The root overlay adds
zero molecular calls; five real-ledger tests pass, with no skipped checks.
The prior archive replay has six passing tests. An
[independent review](INDEPENDENT_REVIEW.md) found no blocking issue.

The smaller pool needs12 rather than20 solvent cells when all three geometries
are distinct. This is40% fewer nominal solvent calls, not a measured40% latency
gain. One recovered Ca proposal reaches the existing movement boundary. Native
solver qualification limits remain; the failed separate solvent-force experiment
does not enable force-driven composite optimization here.

This supports an opt-in minimal adaptive candidate for improved individual-fold
robustness, with a remaining wrong fold, an inconclusive fold and unavailable
cases. Production defaults and explicit DFT remain unchanged. No broad affinity
validation, independent biological test set, or new physiological assignment is
claimed.

## Reproduce the completed overlay

Run from the repository root; use a new output filename because results are
immutable. This reads existing artifacts and performs no molecular calculations.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 scripts/adaptive_minimal_recovery_compare.py \
 --base workspaces/adaptive_minimal_pool_20260923/COMPARISON_v1.json \
 --recovery workspaces/adaptive_origin_recovery_20260923/run_v2/pool/final_collection.json \
 --output workspaces/adaptive_minimal_pool_20260923/COMPARISON_recovered_replay.json
```

Authoritative completed output:
`workspaces/adaptive_minimal_pool_20260923/COMPARISON_recovered_v1.json`, SHA256
`47c30d3f4b50a0fa9da1db751900d3604d57b1648201214c404be7a2b361dafd`.
The exact reference and source pins are in the independent review. The earlier
archive-only report remains unchanged.
