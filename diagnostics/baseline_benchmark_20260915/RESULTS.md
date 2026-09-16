# Baseline benchmarking results

## Outcome

**Retain the baseline as the reference/default.** All 26 endpoint calculations completed, with no retries or failures. No production rescore, threshold fit, or default change occurred.

- **Frozen GGR original v2:** S = -1.183146 kcal/mol, passes its preregistered Ca-direction test.
- **Repaired GGR v3:** S = +3.057473 on the same ranking gauge; the sign reverses. The chemically justified amide repair does not show predictive improvement on this control. V3 has no inherited zero threshold or calibrated decision bands; this is a development observation, not a second confirmatory pass/fail result.
- **Frozen 1KB0 PQQ:** S = 13.254600, below the released Ca upper band 14.857129; class-transfer test passes. This is not a direct affinity measurement.
- All six repairs shift R by +2.815552 to +4.240618 kcal/mol, preserving formal charges and source donor geometry. Repair/original endpoint-pair assigned rank-time ratios are 1.0561–1.1378 on the same node. Affordability is promising; improved discrimination is not established.

**Cost:** Slurm job 1199299, 807 s (13m27s), 64 allocated CPUs, **51,648 allocated core-s**, recorded CPU time 40,778 s, peak batch RSS 8,280,788 KiB (7.897 GiB), zero GPU use. Workflow timing excludes about two seconds of batch startup/finalization and is recorded separately below; do not add these overlapping measurements. Endpoint timing includes launch/startup. No repeat calculations of the existing 27 calibration/holdout cases or four aquo endpoints were needed.

**Checks:** 23 existing real-artifact regression/preparation tests plus four new benchmark/reference/receipt tests passed (27 total, no skips). All 26 new ORCA calculations converged and terminated normally; input, coordinate, output, runtime and execution-receipt hashes were checked. These execution checks do not validate prediction accuracy.

Canonical PQQ remains useful, while its perfect motif/composition separation limits claims of additional DFT information. This expansion contains only one direct same-assay direction control. Broader composition-challenging direct La/Ca labels remain necessary; no new cases were invented or relabeled. Aequorin and parvalbumin remain separate evidence strata and ordered vectors.

The global route is [archived](../global_representation_20260915/ARCHIVE.md). Neither global nor APBS environmental scores enter this benchmark. No response correction was enabled. The v3 chemistry repair was already implemented before this phase; this phase executed it, added comparison/provenance/reporting, and measured its effect. No new scientific protocol ID was created.

Commands: [COMMANDS.md](COMMANDS.md). Unrounded compact ledger: [RESULT.json](RESULT.json). Existing inputs and historical results remain untouched.

Status: complete; 26/26 endpoints complete.

S uses the verified symmetric CN8 ranking gauge. Larger values are more La-like on each protocol scale. Generic and repaired protocols have no calibrated classification bands. Only exact fixed-core PQQ inherits its released bands.

| Case | Representation | S (kcal/mol) | Decision / frozen test |
|---|---|---:|---|
| aequorin_1sl8_EF1 | baseline | -12.584212 | uncalibrated_protocol |
| aequorin_1sl8_EF1 | repaired | -9.219111 | uncalibrated_protocol |
| aequorin_1sl8_EF3 | baseline | -1.998532 | uncalibrated_protocol |
| aequorin_1sl8_EF3 | repaired | 2.241352 | uncalibrated_protocol |
| aequorin_1sl8_EF4 | baseline | 5.273524 | uncalibrated_protocol |
| aequorin_1sl8_EF4 | repaired | 8.569449 | uncalibrated_protocol |
| carp_parvalbumin_4cpv_CD | baseline | -11.590216 | uncalibrated_protocol |
| carp_parvalbumin_4cpv_CD | repaired | -8.774664 | uncalibrated_protocol |
| carp_parvalbumin_4cpv_EF | baseline | 14.419733 | uncalibrated_protocol |
| carp_parvalbumin_4cpv_EF | repaired | 18.334090 | uncalibrated_protocol |
| ggr_1glg_GGR | baseline | -1.183146 | pass |
| ggr_1glg_GGR | repaired | 3.057473 | uncalibrated_protocol |
| pqq_1kb0 | fixed_core | 13.254600 | Ca-supported; pass |

## Chemistry repair effect

Differences below are repaired minus original R; the aquo offset cancels.

- aequorin_1sl8_EF1: 3.365100 kcal/mol
- aequorin_1sl8_EF3: 4.239884 kcal/mol
- aequorin_1sl8_EF4: 3.295925 kcal/mol
- carp_parvalbumin_4cpv_CD: 2.815552 kcal/mol
- carp_parvalbumin_4cpv_EF: 3.914357 kcal/mol
- ggr_1glg_GGR: 4.240618 kcal/mol

## Interpretation limits

GGR supplies one direct same-assay Ca-favoring control. Its frozen sign criterion applies to the original v2 record; repaired v3 is a development comparison. Aequorin EF1/EF3/EF4 is an ordered vector against protein-level evidence, not three labeled independent sites. Parvalbumin CD/EF is supporting cross-study evidence. 1KB0 tests PQQ class transfer, not measured relative affinity. Structural sites are grouped by protein. No flexible classifier or threshold was fitted. All unavailable endpoints remain in the denominator.

Typed donor composition, denticity, core charges and mean donor distance are in the JSON ledger. These are descriptive comparators without a fitted decision rule. Canonical calibration already separates by composition; this panel cannot alone prove incremental DFT value.

Baseline defaults and historical results are unchanged. The APBS and global challengers are excluded. No environmental or response correction is supplied.

## Measured execution cost

- Job 1199299: 804.859 workflow seconds, 64 allocated CPUs, 51510.959 allocated core-seconds. Scheduler totals and peak memory are recorded separately after termination.

Unrounded energies and pinned receipts: `/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/baseline_benchmark_20260915/run_v1/collection_1199299.json`.
