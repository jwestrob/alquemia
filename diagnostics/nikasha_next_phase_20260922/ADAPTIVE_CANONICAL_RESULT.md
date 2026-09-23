# Adaptive completion: full canonical calibration is now available

The same four force-selected angular motions now produce valid candidates for
all30 original sources and both metals. One fixed optimizer unit scaling resolved
the earlier MMOL1770 failure without changing chemical states or final geometry
limits. Protocol: `common_four_angular_native_OMOL_SLSQP_hartree_units_v2`.

The common-pool scorer cross-evaluated every admitted candidate for both metals.
All30 pools are available:60 fresh cross-MACE evaluations and240 native GFN2
singlepoints completed, with60 candidate MACE evaluations and all existing
three-candidate cells reused. No DFT, refolding or new chemical state was used.

## Fair calibration and remaining question

| Population | Old static bands | New canonical-only bands |
|---|---|---|
| Original25 canonical references |23 correct,1 wrong,1 inconclusive|25 correct|
| Three consumed crystal controls |3 correct|3 correct|

The two unknown PLM examples do not supply labels or set thresholds. Both
mathematical and operational row selection yield the same frozen bands:

- Ca-supported: `R <= -405459.0155077257` model kcal/mol.
- La-supported: `R >= -405456.46881754015` model kcal/mol.
- Between them: inconclusive. Gap:2.5466901855543256 model kcal/mol.

`R = E_Ca - E_La`; this model-dependent electronic contrast has no compatible
aquo reference and is not a binding free energy. Calibration uses exactly the
released25 members and unchanged extrema/minimum-gap rule. The reference was
frozen before new primary225 search execution. This establishes retained
reference fidelity, not improved fold robustness or independent validation.

The Q9Z4J7 old-band wrong call and Q88JH5 abstention disappear under the new
canonical reference. This does not establish improvement: those are calibration
members. Full225 transfer, using this frozen reference, is the next comparison.

## Physical and cost limits

All60 optimizations succeeded with nonpositive native MACE work.27 endpoints
reach an angular/displacement boundary; these are finite constrained proposals,
not unconstrained stationary minima. The common composite pool selects58
adaptive candidates among60 metal endpoints; candidate frequencies are not
populations. Individual Ca/La works and all geometry selections remain in the
machine-readable inspection.

Proposal workers performed1,475 fresh native MACE energy/force evaluations;
the60 cross-scoring calls bring this stage to1,535 fresh MACE evaluations.
All60 exact origin evaluations were reused. Each count comes from the completed
worker summaries, including trial evaluations rather than counting only final
optimized endpoints.

Actual Slurm allocation for proposal generation plus common-pool evaluation:
60,864 CPU-seconds and548 GPU-seconds. This includes both proposal jobs and all
five scoring jobs, including orchestration/collection within allocations; it
excludes local preparation and reused historical work. These are development
batch costs on the two recorded CPU hosts, not a newly measured single-source
production latency. No production promotion occurred.

## Artifacts

- `workspaces/adaptive_completion_20260922/original30_v1/TERMINAL_SUMMARY.json`
- `workspaces/adaptive_completion_20260922/original30_pool_v1/final_collection.json`
- `workspaces/adaptive_completion_20260922/original30_pool_v1/REFERENCE.json`
- `workspaces/adaptive_completion_20260922/original30_pool_v1/inspection.json`
- `workspaces/adaptive_completion_20260922/original30_pool_v1/ACCOUNTING.txt`

Jobs1209963,1209968,1209983,1209985,1209987,1209989,1209991 are complete.
The prepared primary225 manifest retains all225 sources/450 endpoint statuses;
20 inherited unavailable pools remain explicit. Execution staging/parallelism
is recorded separately by the adaptive agent.

The unchanged410 eligible searches were submitted in four disjoint shards as
1210021–1210024 under `primary225_v2_sharded/manifest.json`. This changes execution
parallelism only; the older unlaunched manifest is preserved. No transfer result
was used to choose this calibration or any scientific parameter.
