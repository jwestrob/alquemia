# Khoury author-domain benchmark: two supporting successes, RTX failure

2026-09-18. The frozen existing masked MACE descriptor passes **6/9 declared
domain-mean comparisons**, covering **two of three new domain groups** against
the three consumed GGR structures. RTX fails all three comparisons. These are
qualified protein-level/cross-readout comparisons, not nine independent
biological tests or fitted La/Ca affinity predictions. Production is unchanged.

| Author domain | Sites | Declared mean, model kcal | Mean−GGR1GLG | Mean−GGR2FW0 | Mean−GGR2FVY | Passes |
|---|---:|---:|---:|---:|---:|---:|
| A0A7 | 6 | 49.234436 | +25.247014 | +4.350308 | +3.258935 | 3/3 |
| HEW5 | 8 | 70.120351 | +46.132929 | +25.236223 | +24.144850 | 3/3 |
| RTX | 8 | 2.777993 | −21.209429 | −42.106134 | −43.197508 | 0/3 |

The declared screen requires a margin greater than 0.02 model kcal. The mean
uses every modeled site with equal weight. No threshold, site, water inventory,
charge state or aggregation rule changed after scoring. No PQQ band or universal
zero was applied. All 22 sites and 44 endpoint forwards completed; none is omitted.

## All ordered site scores

Site order is the original author metal-chain order. Each state substitutes
only the selected ion; all other ions in that domain remain Ca2+.

| Domain | Ordered site vector, model kcal |
|---|---|
| A0A7 B–G | 61.401780,56.767042,30.462803,48.987096,53.803890,43.984004 |
| HEW5 B–I | 60.969255,74.129786,74.606595,66.088231,73.017751,69.436125,66.968435,75.746626 |
| RTX B–I | 4.327496,5.294125,15.785283,−2.885218,−7.706783,−10.360517,6.714642,11.054914 |

Every RTX site is below every GGR comparator. This is a failure under the
declared qualified comparison; it is not explained away by selecting another
RTX site. A0A7's favorable mean also does not imply every A0A7 site individually
exceeds every GGR structure. Site-level experimental La/Ca labels are absent.

## Evidence and representation limits

SI Table S6 of [Khoury et al.2025](https://doi.org/10.1039/D5SC02315G) reports
La ITC Kd17±2µM(A0A7),5.2±1.7µM(HEW5),40±4µM(RTX),
with 95% confidence intervals. Calcium≈2000,750,250µM are CD concentrations
inducing folding, not fitted Kd. ITC and CD differ in buffer, salt and protein
concentration; the reported multisite affinities cannot label individual sites.
Thus these results support an expanded exploratory benchmark, not a new
direct same-assay gold-standard panel. No within-lanthanide or physiological
occupancy claim follows.

All three selected domain sequences match SI Table S3 exactly. Experimental
cleaved constructs retain an N-terminal MPVP scar absent from the author models.
No scar was modeled or silently omitted from the description. The supplied
Ca-conditioned coordinates were preserved; these are not measured La geometries.
The paper's Fig. 1 caption mentions RTX reference 5CVW while describing the set as
AF3 predictions; this pilot uses the actual author-supplied s008 coordinates,
not downloaded 5CVW coordinates or a claim that s008 is experimental geometry.

All seven author PDBs are inventoried. Their protein lengths/Ca counts are:
s00187/6,s00273/2,s003117/8,s00443/1,s005151/10,s00664/4,s008152/8.
The four unscored models remain inventory-only; weak XO/FRET responses were
not converted into calcium-preference labels. Original heavy atoms map exactly
through metal-chain normalization. Native pH6 hydrogen preparation adds no
water or missing heavy atoms; all three sources already contain OXT. HEW5 has
one HIP and RTX five HIP under the fixed native rule. Ca/La share each recorded
protonation/H geometry. No score-conditioned preparation rescue occurred.
The retained GGR2FW0/2FVY comparators use archived pH7 protonation, whereas
these new domains use pH6. Thus this pilot also retains that preparation
difference; no matched-pH control was run, and no unique cause is assigned
to RTX's failure. This does not change the predeclared result.

## Validation, cost and reproducibility

Three real-fixture tests passed in 25.400 s, no skips: unchanged legacy parvalbumin
replay, exact author-atom mapping/corrupted-copy rejection, and all 22 paired
states/background inventories. All 44 actual endpoint accounting checks pass.
The repeated all-Ca state agrees exactly within each domain despite reordered
atoms (measured spread 0.0; frozen numerical tolerance 0.01 model kcal).

Job 1201351 completed on one A5000 with 16 CPUs and 64474 MiB host allocation:
**799 wall/GPU-allocation seconds;12784 allocated core-seconds;852.837 reported
CPU-seconds**. The44 native forwards took 198.585446 s in total; peak GPU allocation
was 3243407360 bytes and peak worker host RSS 1421544 KiB. Initial source/preparation
generation recorded 517.424638 wall seconds and 449.499209CPU seconds. Two transient
preflight failures caused by editing the plan during preparation are retained;
the original plan was restored and the same inputs recovered before inference.
No extra model/DFT/folding/relaxation calculation was needed. Recovery housekeeping
was not separately timed. The read-only report took 549.30 wall seconds and
476.58 CPU seconds, with no new inference; its resource receipt is retained.

New preparation policy:
`omol_author_AF3_domain_multisite_pH6_fixed_Ca_ff19sb_v1`.
The existing descriptor/checkpoint/factorization and old multisite default are
unchanged. No baseline result, default, calibration or historical record changed.

Products: `workspaces/mace_omol_20260917/khoury_author_domains_v1` (sources and
original attempts),`khoury_author_domains_recovery_v1/manifest.json` (complete
22-site inventory;SHA adcc267275f599dfff1097a5cac2af21bc1ed77dfb64a6406e1c94badf5dbd84),
and`khoury_author_domains_report_v1/result.json` (full unrounded comparison).
The recovery launcher directory holds actual submission/accounting receipts;
`khoury_author_domains_recovery_v1/cost_with_report_v1.json` retains exact measured costs.
each site retains its frozen worker, native readouts and attempts. Exact commands
are in KHOURY_BENCHMARK_COMMANDS.md; selection/aggregation were declared in
KHOURY_BENCHMARK_PLAN.md before any score. No additional pilot is launched here.

**Recommendation:** add these three qualified domain groups to the research
ledger and retain the production baseline. The new evidence broadens what the
existing MACE descriptor handles, while RTX supplies a useful independent
development failure. It does not establish a broadly reliable replacement or
validate the unchanged DFT baseline on these new domains.

## Figure export

A plot of all22 ordered site scores, the three fixed domain means, and all three
GGR structures is available as [PDF](../../workspaces/mace_omol_20260917/khoury_benchmark_figure_v1/benchmark.pdf),
[SVG](../../workspaces/mace_omol_20260917/khoury_benchmark_figure_v1/benchmark.svg), and
[PNG](../../workspaces/mace_omol_20260917/khoury_benchmark_figure_v1/benchmark.png).
The shaded interval is the observed GGR structure range, not a classification
band. The plotting script and source hash sit beside the exports; no new
selection, model calculation or statistical comparison was added.
