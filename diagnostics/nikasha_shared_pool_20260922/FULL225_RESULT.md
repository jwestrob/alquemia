# Full225 shared-pool result

**Retain the released static context composite. Giving both metals the same finite geometry pool did not improve discriminatory performance or earn routine deployment.** All 406 new MACE and 812 new GFN2 evaluations completed successfully. The 20 unavailable sources were inherited; this run added no new failures.

## Actual utility

Counts are **correct / wrong / inconclusive / unavailable**. All 225 sources and 100 overlapping La-conditioned triples remain in the denominators. Mathematical and operational pool variants give identical decisions throughout this challenge; their raw values remain separate.

| Method | All225 | La100 | Ca125 | Strict100 triples |
|---|---|---|---|---|
| Preserved DFT | 199/3/6/17 | 96/1/1/2 | 103/2/5/15 | 94/0/0/6 |
| Native MACE core | 199/0/9/17 | 95/0/3/2 | 104/0/6/15 | 94/0/0/6 |
| Released static context composite | 203/2/2/18 | 97/0/1/2 | 106/2/1/16 | 94/0/0/6 |
| Own proposal, old bands | 201/1/3/20 | 97/0/1/2 | 104/1/2/18 | 94/0/0/6 |
| Own proposal, its calibration | 197/2/6/20 | 95/0/3/2 | 102/2/3/18 | 92/0/2/6 |
| Shared pool, old bands | 201/1/3/20 | 97/0/1/2 | 104/1/2/18 | 94/0/0/6 |
| Shared pool, its calibration | 198/2/5/20 | 95/0/3/2 | 103/2/2/18 | 92/0/2/6 |

On the **same 205 available sources**, static context gives **201/2/2**, while the newly calibrated pool gives **198/2/5**. Three correct calls become inconclusive; no wrong call is recovered. Two additional static-correct sources remain unavailable to this entire proposal branch because its earlier proposal solvent calculations failed.

The three matched transitions are C5AXV8 La-conditioned samples 2 and 3, and Q88JH0 Ca-conditioned sample 0. The same two C5AXV8 three-fold subsets become inconclusive: sample sets {0,2,3} and {2,3,4}. The released context, preserved DFT and native core retain all 94 available triples correctly; the new pool retains 92, with two inconclusive and the same six unavailable. These 100 overlapping subsets repeat 25 biological groups.

The old-band pool exactly reproduces the earlier own-proposal decisions: it defers one static wrong call (A0A3F2YLY8 Ca sample 1) to inconclusive, without increasing matched correct calls. This is an unchanged old-band transfer check, not the pool’s own calibrated classifier.

Compared with the previous separately calibrated own-proposal model, the pool gains one correct call (A0A3F2YLY8 Ca sample 0), **with exactly unchanged raw R**. That gain comes from the distinct canonical calibration boundary; it is not evidence of a better energy for that fold.

## Frozen reference and strict aggregation

- The pool reference was frozen before full-fold execution: `remaining26_v1/REFERENCE.json`, SHA256 `f96b2425202c0c24f11312f28109286cb88b5fb520a75a16d2a9e4ee9dbd0911`.
- Both variant-specific bands are Ca_max −405462.0237183686 and La_min −405456.53376105055 model kcal/mol. Only the original 25 canonical members set these extrema; no fold result changed the bands.
- Strict La4: static, DFT and both pool variants each 23 correct / 2 unavailable. Strict Ca5: static 21 correct / 4 unavailable, pool 19 / 6; preserved DFT 21 correct / 1 inconclusive / 3 unavailable. Equal-arm median averages: static 20 correct / 5 unavailable, pool 18 / 7, DFT 21 / 4. Missing members are never dropped.
- On the same complete groups, operational La4 range medians change 9.08319→10.29877 model kcal/mol (13 groups decrease, 10 increase); Ca5 changes 12.12867→10.05755 (10 decrease, 9 increase). This mixed structural spread is not an accuracy gain or thermal uncertainty estimate.
- Released context has the most correct individual-fold calls here. Native core gives fewer correct calls but zero wrong calls; this report does not claim a universal winner for every error tradeoff.

## Execution and cost

- GPU 1209870: 406/406 MACE calls. CPU 1209871/1209872/1209874/1209880: 812/812 GFN2 calls. All terminal successful; no retry or new scientific computation during reporting.
- Pending shard-2 job 1209873 was canceled before any attempt and replaced by the same finite manifest on a second CPU host as 1209880. No task duplication; cancellation and replacement receipts retained.
- Incremental allocation: **169,984 core-seconds and 116 GPU-seconds**. Source archive preparation separately recorded 63.650070 wall seconds. This excludes prior origin/proposal work and canonical calibration. The mixed CPU hosts do not support a matched-hardware speed claim.
- Inherited unavailable cases: 17 unsupported source preparations; one original A8R3S4 Ca sample 3 solvent endpoint; earlier proposal endpoints for P12293 Ca sample 4 and Q60AR6 Ca sample 1. No substitution with static scores.

## Artifacts and recommendation

- Final actual collection: `workspaces/nikasha_shared_pool_20260922/primary225_v1/final_collection.json`.
- Full matched methods, every source, strict pool and triple: `workspaces/nikasha_shared_pool_20260922/primary225_v1/final_comparison.json` and `.md`.
- Compact transitions, costs, spread and raw shifts: `workspaces/nikasha_shared_pool_20260922/primary225_v1/final_summary.json`; scheduler record `scheduler_accounting_final.txt`.
- Existing comparator validated saved matrix selection/algebra and original source membership. Reporting checked exact old-band decision parity with own proposals, mathematical/operational decision parity and all 225/75/100 denominators.

**Recommendation:** keep the released static method for routine use; do not promote this restricted shared-pool branch. The common-pool implementation closes the fairness question for these candidates and can support separately declared adaptive-candidate tests. It has not demonstrated an accuracy improvement itself. Baseline and defaults remain unchanged.
