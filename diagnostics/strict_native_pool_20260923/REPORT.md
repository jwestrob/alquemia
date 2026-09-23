# Uniform strict native scalar test: qualified on32 consumed sources

Fresh strict starts and strict continuations agree across **384/384 exact cells**, all96 same-geometry contrasts and all32 complete three-geometry pools. Both branches give **25/25 canonical,3/3 crystal and4/4 noncanonical calls**, with either the earlier bands or their own canonical-only references. The earlier loose-cold panel gave31correct/1inconclusive. Production/defaults remain unchanged.

## Numerical finding

| Frozen check | Maximum observed difference | Allowed |
|---|---:|---:|
| Per-cell energy, seeded minus fresh | 0.0000103489 kcal/mol | 0.1 |
| Same-geometry Ca−La composite contrast | 0.0000100496 kcal/mol | 0.2 |
| Three-geometry pooled contrast | 0.00000282377 kcal/mol | 0.2 |

All384 fresh calculations confirmed SAD with explicit NoAutostart and no starting seeds. All384 seeded results (374new,10exactreused) confirmed XTBRESTART from the exact original loose-cold pair. Every accepted output verifies native mixer, effective TolE1e-10Eh, state, parameter export and charge sanity. No failed calls or retries occurred. The supplementary pass2-seeded strict20 records were not substituted into either branch.

Printed RMS-density limits pass383/384fresh and384/384seeded; MAX-density passes320/384 and327/384. Native energy stopping and scalar repeatability are supported here. Analytic-force consistency, electronic ground-state uniqueness, equilibrium populations and broad biological accuracy are not established. Do not convert these tests into a force qualification.

## Discrimination and interpretation

Only one of64 endpoint candidate choices changes relative to loose cold: La for Q4W6G0 Ca-conditioned sample2 selects adaptive_La instead of adaptive_Ca, identically in both strict branches. Its raw R changes from−405463.62238655 to−405471.60959642 modelkcal/mol and becomes Ca-supported under both old and frozen new bands. The largest corrected scalar component is its La/adaptive_Ca/vacuum cell:−11.15558121kcal/mol relative to loose cold. This is the geometry-sensitive numerical anomaly identified earlier; no favorable geometry or initialization was selected to repair it.

The other three noncanonical sources retain their calls: P38539La2 Ca-supported; A0A3F2YLY8Ca1 and A0ACD6B9F2Ca4 La-supported. The crystals1H4I/4MAE/1KB0 retain Ca/La/Ca. These are consumed, correlated development sources, not32 independent blind biological observations.

Fresh canonical-only bands are Ca_max=−405463.71328194055 and La_min=−405456.37993968255; gap**7.3333422580**modelkcal/mol. Seeded bands remain separately recorded (gap7.3333408788). The prior gap was7.3253888497; the small spread change is not a substantive chemistry improvement. Mathematical and operational selections coincide here. Exact unrounded scores and both references are retained; no threshold fitted to the four folds/crystals.

`REFERENCES.json` SHA256:`171bd31d466ff97ef6073ce23286af8519ef23690f7d6ad23bb444cab30ecd5c`. Fresh branch schema:`branches.fresh.variants.{mathematical,operational}.bands`.

## Actual cost

Job**1211270** completed758newscalarcalls plus10explicitreuses in**510s on32CPUs =16320allocatedcore-s**,64GiB requested,zeroGPU. This is development validation of both starts.

- Fresh384call executor:226.770720s;7256.663033allocatedcore-s during that interval; summedORCAruntimes6407.747s.
- Seeded374newcall executor:156.152765s;4996.888470allocatedcore-s; summedORCAruntimes4469.204s. Its historical loose-cold seed cost is additional.
- Remaining job time covers validation, copying, collection and other wrapper overhead. This overhead is included in total allocation. Local preparation/tests/report are outside that allocation.

The fresh scalar stage needs one call percell. The seeded route requires a loose cold call plus one strict continuation. No MACE/geometry-search costs were incurred here; this is not the full scanner runtime. Four real-artifact tests passed in16.253s,zero skips.

## Next decision

The fresh strict protocol earns a fixed-geometry structural-transfer test under its frozen reference. It repairs a demonstrated numerical defect without changing chemistry or adding a routine second electronic calculation. The separately authorized four-source comparator and full225transfer retain the original ledgers and will determine broader utility. No additional starts, stronger thresholds or force model are introduced by this result.

Exact results/costs:`workspaces/strict_native_pool_20260923/run_v1/{COMPARISON,REFERENCES,COSTS}.json`. Recipe/qualification/command pins are in `ARTIFACTS.json` and `COMMANDS.md`.
