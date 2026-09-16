# Baseline readiness for PLM PQQ candidates

2026-09-16. Read-only review of existing results, requested before expanding use.
No folds, protonation retries, energy calculations or threshold changes.

## What remains successful

The exact `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3` calibration
separates **25/25 controls (11 La-class, 14 Ca-class)** with an 8.602932 kcal/mol
gap. Frozen crystal transfer passes for 1H4I (S 7.517844) and 4MAE (38.089733).
The later 1KB0 Ca-class test passes (13.254600). These are separate evidence
sets, not 28 independent biological validations. Motif, charge and composition
already separate the original calibration; it does not isolate DFT's added value.

The authoritative supported bands are Ca S <= 14.857129202922806 and La
S >= 23.460061205609236 kcal/mol, with the open interval indeterminate. These
are empirical supported regions, not probabilities or a measured dual-affinity
range. Keep the earlier midpoint criterion distinct from the released bands;
do not switch decision rules to force a PLM call.

## PLM ADH9 results actually available

| Gene | Protenix S | Reviewed AF3 S | Interpretation |
|---|---:|---:|---|
| 32301_3 | 15.503190 | 19.984118 | Both indeterminate; source-model shift +4.480928 |
| 4380_6 | 21.450083 | unavailable | Protenix indeterminate; AF3 preparation invalid |

Units: kcal/mol on the released reporting gauge. The two Protenix pairs use the
correct unchanged fixed-core helpers, native ORCA 6.1.1 r2SCAN-3c/CPCM(Water),
complete oxidized PQQ, 78-atom cores, La/Ca charges -2/-3 and singlets. Existing
receipt/geometry reviews accept them. The AF3 4380_6 source protonation produced
13 H–H overlaps below 0.5 Å (minimum 0.058626 Å); its raw score is withheld.

The original selected panel has **six proteins: two valid pairs, four unsupported
preparations**. Three fail the frozen Asp–cation contact gate, one the CN gate.
Do not report all six as inconclusive scored predictions. This review concerns
that identified ADH9 panel, not every historical PLM application.

## A specific calibration coverage gap

Both scored PLM sites have the additional acidic D+2 residue and a **Lys** partner
to the catalytic Asp. The frozen `core_map.tsv` shows:

| Calibration class | Arg partner | Lys partner |
|---|---:|---:|
| La | 11 | 0 |
| Ca | 12 | 2 |

The Lys Ca controls are Q4W6G0 and Q8GR64. The 4MAE La-class crystal also has
Arg. The later 1KB0 control has Lys but is Ca-class without the extra acidic
residue. Thus the reviewed positive controls do not cover the PLM combination
of extra Asp plus Lys. The protocol permits that combination; successful input
preparation does not establish calibrated transfer to it.

This is a **coverage observation, not evidence that Lys causes a score shift**.
Do not mutate a residue, remove it, change a gate or adjust a threshold to obtain
a desired call. XoxF-like annotation alone is not an experimental La/Ca label.

## Practical judgment and next scientific question

Retain the baseline for qualified PQQ screening with explicit indeterminate and
unsupported outcomes. Existing evidence supports the canonical controls; it does
not establish reliable decisive calls across divergent PLM PQQ proteins, or broad
La/Ca discrimination outside that family. The known generic GGR/alpha limitations
remain separate and have not been solved by the MACE proposal.

The most targeted missing benchmark is an independently characterized La-using
PQQ protein with the **extra Asp + Lys partner** architecture, prepared by the
same protocol. Finding/checking that control is a proposed next analysis; no new
literature screen, structure prediction or calculation was run in this review.
A controlled Arg/Lys perturbation could later examine mechanism, but would be a
development experiment, not experimental validation of the PLM proteins.

## Source records

- [Calibration](../pqq_pmdh_fixed_core_calibration_20260914/RESULT.md),
  [machine-readable authority](../pqq_pmdh_fixed_core_calibration_20260914/result.json),
  [core map](../pqq_pmdh_fixed_core_calibration_20260914/core_map.tsv).
- [Frozen crystal transfer](../pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/result/HOLDOUT_RESULT.md).
- [1KB0 and generic benchmark](../baseline_benchmark_20260915/RESULTS.md).
- [PLM fixed-core preparation/results](../plm_adh9_fixed_core_20260916/README.md).
- [Reviewed AF3 comparison](../plm_adh9_af3_20260916/README.md).
- Native PLM records under
  `/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/adh9/`:
  `execution/results.json`, `execution/terminal_review.json`,
  `af3_comparison/reviewed_score_comparison.tsv`, and
  `af3_comparison/terminal_review.json` superseding the raw AF3 PASS.

See `BASELINE_CHECKS.json` for this review's artifact hashes, archived parsing,
classification replay and descriptor counts. No new accuracy estimate is fitted.
