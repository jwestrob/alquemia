# Lanmodulin ortholog REE profiles: unsealing record and first-look plan (frozen 2026-11-10)

Lead agent: Claude Fable. Written before any outcome value was read. Only the header row and the ReadMe sheet of
the workbook have been read so far (column names; no values).

## Unsealing record

- **Authorized by Jacob on 10 November 2026:** "Unseal the data as you like; just mark that you did it."
- **Source:** `/groups/banfield/projects/environmental/sr/srvp2020/Jacob/spicy_lams/41589_2026_2176_MOESM3_ESM.xlsx`
  (Supplementary Data 1 of Nat. Chem. Biol. s41589-026-02176-3; sha256 `39193323…99f1` per the 1 August
  selection record), sheet `Data`, 616 ortholog rows.
- **Unsealed now, for the 600 rows outside the reserved panel:** raw desorption by ICP-MS (`T:BL`, 15 elements ×
  3 replicates) and normalized logD mean and SD (`BT:CW`).
- **Stays sealed:**
  1. Every outcome cell (`T:BL`, `BT:CW`, `CX:CZ`) of the 16 reserved proteins of
     `spicy_lams_8family_double_difference_sentinel_v1` (8 families, pairs p01–p08). The extraction script
     drops those rows by accession before anything is stored or printed.
  2. The published cluster and UMAP columns (`CX:CZ`) for all rows. Unneeded here.
  3. The article's results sections and the derived preprint. Unread, because they name individual orthologs.
- **Status of the 1 August sentinel (v1):** unchanged. Its fixed-cluster energetic predictions were never
  locked, and its 16 rows remain sealed. A sequence-based protocol (v2) will be written and hash-locked from
  the 600 unsealed rows before the 16 rows are read, and scored with the v1 statistics.
- **Marks:** this file; `diagnostics/spicy_lams_unseal_20261110/UNSEAL.md` in the repository; a notice file in
  `spicy_lams_matched_homolog/`; a dated block in the repository's `AGENTS.md`; the vault `CONTEXT.md` and a
  capture.

## Why open it

Every modelling question of the past week ended at a benchmark too small to answer it. This table is the only
protein-scale measurement set in the project (616 proteins × 15 elements × 3 replicates). Jacob needs to know
whether asking a laboratory to measure further proteins could be worth it. That depends on three things this
table can show: whether natural lanmodulins differ reproducibly, by how much, and whether sequence predicts the
differences.

## Data handling (fixed)

- Assay (methods already on file): SpyTag–lanmodulin on SpyCatcher glass, 100 ppb of each rare earth at pH 3,
  15 min, wash, acid desorption, ICP-MS. Total captured metal is a few percent of the load, so the captured
  mixture reflects competition among elements.
- Replicate profile: `l[p,i,k] = log10(desorbed ppb)` for protein p, element i, replicate k; non-positive or
  missing values are missing. Centred replicate profile: `l` minus its mean over the 15 elements. Mean profile:
  the mean of the centred replicate profiles. Energies: −1.3643 × log10 units (kcal/mol, 298 K).
- **QC set:** total-metal mean ≥ 20 ppb and RSD ≤ 50% (the blind rule of 1 August) and all 15 elements positive
  in at least two replicates. All analyses use the QC set unless stated.
- Sanity check: correlation of my centred mean profile with the workbook's normalized logD mean (centred).
- Size axis: Shannon nine-coordinate radii as in the 10 November gap analysis; Y at 1.075 A.

## Step 0: the frozen predictions of this morning (scored first)

`FROZEN_PREDICTION.md` of the gap analysis (sha256 `3f9b2507…13ec` for the blind text):

- **P1:** pooled gap between Gd and Tb below +15 mA (cubic polynomial per ortholog, one shared gap).
  Contradicted by a 95% interval above +18.5 mA. Clearly negative (below −10 mA) reads "chelator-referenced".
- **P2:** yttrium minus the ortholog's cubic polynomial at 1.075 A averages −0.05 to +0.25 kcal/mol and Dy or Ho
  is the closest lanthanide. Revised expectation after the public lanmodulin series (not blind): +0.4 to +0.9,
  closest to Er or Tm.

## Step 1: first-look readings (frozen)

Per ortholog and per replicate, fit the 14 lanthanides with orthogonal polynomial components in radius
(constant, linear, quadratic, cubic; Legendre-type on the La–Lu range) and record the yttrium offset (Y minus the
cubic at 1.075 A).

- **A1. Reproducible variation.** One-way intraclass correlation over the three replicates, across orthologs,
  for the linear, quadratic and cubic components and the yttrium offset. "Reproducible variation exists" if the
  linear and quadratic components both have ICC ≥ 0.5 with a bootstrap 95% lower bound ≥ 0.4 (bootstrap over
  orthologs).
- **A2. Smoothness.** Residual mean square of the cubic fit to each mean profile against the replicate noise of
  that mean (pooled). "Profiles are smooth" if the ratio is at most 2.
- **A3. Size of natural differences (descriptive).** Noise-corrected SD across orthologs and the 5th–95th
  percentile range, in log10 units, for Lu−La, Dy−Nd, Nd−Pr, Dy−Tb, Tb−Gd, Y−Dy and Y−Ho.
- **A4. Does sequence predict profiles in unseen families? (primary)** Leave-one-family-out over NCBI families
  with at least three QC orthologs; all other orthologs always train. Target: the mean centred 15-element
  profile. Predictors:
  - B0: training mean profile.
  - B1 (primary): profile of the training ortholog with the highest global sequence identity (BLOSUM62, gap
    open −10, extend −0.5; identity over aligned columns including gaps).
  - B2: identity-weighted mean of the five nearest training orthologs.
  - M1: ridge regression on EF-loop features (below), penalty chosen by inner cross-validation.
  - M2: ridge regression on mean-pooled ESM-2 650M embeddings (cached weights, run locally).
  Scores, families weighted equally: skill = 1 − MSE(model)/MSE(B0), and the mean cosine between (prediction −
  training mean) and (observed − training mean). Bootstrap over families.
  "Sequence predicts selectivity across families" if B1 has skill above 0 with the interval excluding 0 and
  mean cosine ≥ 0.3. The other predictors are reported beside it; a better one becomes the candidate for v2.
- **A5. How close must a neighbour be? (descriptive)** Ten-fold random cross-validation of B1; skill and cosine
  by identity to the nearest training ortholog (< 50, 50–70, 70–85, 85–95, > 95%).
- **A6. Which sequence features track which component? (descriptive)** Spearman correlations of the linear,
  quadratic and cubic components and the yttrium offset with: number of EF loops, mean acidic residues per
  loop (loop positions 1, 3, 5, 9, 12), loops with Asp at position 9, loops with Pro at position 2, net charge
  per residue, length. Reported pooled and within family (family means removed), with family-level bootstrap.
- **A7. Capture against loop acidity (descriptive).** Total captured metal against mean acidic residues per
  loop, over all 600 rows. Capture mixes expression with binding, so this is a pointer only.

EF loops: 12-residue windows matching `D.[DNS].[DNSTG]......[DE]`, non-overlapping, kept in order; the count is
compared with the workbook's `EFhands` column, and orthologs where the two disagree are flagged and kept.

## Step 2 (after step 1 is reported)

Lock a sequence model as sentinel v2 and score the 16 reserved proteins once. Apply the locked model to the
SR-VP and Omnitrophota homologs, with identity to the nearest library member beside each prediction.

## What would make further measurements worth asking for (stated now)

1. A1 passes: natural lanmodulins differ reproducibly.
2. A3 shows within-branch contrasts (Nd−Pr, Dy−Tb) or the yttrium offset varying across orthologs by a
   noise-corrected SD of at least 0.1 log10 units: there is something to mine.
3. A4 or A5 shows positive skill at the identity range of Jacob's environmental homologs: there is a way to
   choose them.

If 1 or 2 fails, natural lanmodulin diversity is not worth mining for selectivity. If only 3 fails, sequences
would have to be chosen by diversity alone, at a much lower expected yield.

## Outputs

`scripts/`, `results/` (`run.log`, tables, `summary.json`), `figures/`, `REPORT.md`, `DEVIATIONS.md` if needed.
