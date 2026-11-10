# Addendum to PLAN.md: sample rule and metrics for test 1 (frozen 2026-11-10, before deep-water values were examined)

State of knowledge at freezing: the article text (Shiller et al. 2017) was read; the downloaded spreadsheet was
opened only to see its layout (header and the first twelve data rows, which are surface samples plus two samples
at 210 m and 884 m). No ratio between plume and background had been computed.

## Data

- `data/griidc_download/shiller_ree_griidc.xlsx` (GRIIDC doi 10.7266/N7CC0XQ6): dissolved Y and 14 lanthanides,
  pmol/kg, four cruises.
- `data/griidc_download/shiller_gomri_tracemetals.xlsx` (GRIIDC doi 10.7266/N7MS3QQ5): oxygen, nutrients and
  trace metals for the same samples.

## Readings

- **R1, replication of the article's Figure 2.** May 2010 samples (both cruises) at 800 m or deeper. Plume set:
  the five samples with the lowest La/Yb. Background: every other May 2010 sample at 800 m or deeper. Removal of
  element i: f_i = 1 − mean(plume_i) / mean(background_i); 95% interval by bootstrap over samples in both sets.
  Sensitivity: plume set chosen by lowest La; depth floor 1,000 m.
- **R3, selector independent of the rare earths.** Oxygen deficit of each May 2010 sample at 900–1,400 m =
  baseline(depth) − measured oxygen, with the baseline a straight line in depth fitted to the October 2011
  samples at 800–1,600 m. For each element: slope of ln(concentration) on oxygen deficit; removal per 10 µM
  deficit = 1 − exp(10 × slope); bootstrap interval. This reading guards against the selection effect in R1
  (choosing low La/Yb samples inflates La removal and deflates Yb removal by construction).
- Rate-constant ratios: k_i / k_La = ln(1 − f_i) / ln(1 − f_La).

## Decision terms

- "Nd ≥ La": the point estimate of f_Nd is at least that of f_La in R1.
- "Sm–Gd removal detectable": the mean of f_Sm, f_Eu, f_Gd is at least 0.05 with a bootstrap interval that
  excludes zero in R1, and the corresponding R3 slope is negative with an interval that excludes zero.
- "Vanishes beyond Nd": the mean of f_Sm, f_Eu, f_Gd is below 0.05 or its interval includes zero.
- R1 and R3 must agree for a categorical reading; disagreement is reported as unresolved.

## A physical correction to the plan's lanmodulin prediction (added here, before the data)

PLAN.md predicted the removal pattern of a lanmodulin-like sink from the capture profile alone. That profile was
measured at pH 3 with equal amounts of every element as free ions. In seawater most of each element sits in
carbonate complexes, and the free fraction falls from La to Lu. A site of intrinsic affinity K_i therefore
removes element i at a rate proportional to K_i × (free fraction of i).

- The plan's categorical reading stays as frozen.
- Beside it, a speciation-corrected prediction is reported: k_i / k_La = 10^(P_i − P_La) × α_i / α_La, where P is
  the mean log10 capture profile of the 353 quality-controlled lanmodulins and α is the free-ion fraction from
  published carbonate constants for deep Gulf water (no parameter is fitted to the plume data).
- The same correction applies to any binder, so free-ion fractions are also applied when a natural
  solid-over-seawater ratio is converted to a binder's intrinsic profile in tests 2–4.
