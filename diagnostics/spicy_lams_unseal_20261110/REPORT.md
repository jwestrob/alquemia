# 616 natural lanmodulins: first look after unsealing — 10 November 2026

**Result.** Natural lanmodulins differ reproducibly in how they capture the rare earths, but almost only along
one axis: a tilt between light and heavy elements. Every ortholog separates neighbouring elements, and rejects
yttrium, in nearly the same way. A measured relative at 70% identity or more predicts a protein's profile. On
this evidence, measuring further natural lanmodulin relatives would mostly confirm what the 616-protein library
already shows, and the more distant environmental families lie outside what the library can predict or what this
assay can test.

Jacob authorized the unsealing on 10 November 2026 ("Unseal the data as you like; just mark that you did it").
Record: `diagnostics/spicy_lams_unseal_20261110/UNSEAL.md`. Opened: raw desorption and normalized logD for the
600 workbook rows outside the reserved panel. **Still sealed: all outcomes of the 16 reserved proteins, and the
published cluster and UMAP columns.** `PLAN.md` (sha256 `155dccc5…0363`) was frozen after reading only the
column headers. Lead agent: Claude Fable; Mac, CPU, about five minutes of compute.
Figure: `figures/lanmodulin_ortholog_first_look.png`.

## Data

- SpyCI-LAMBS: immobilized lanmodulin, 100 ppb of each of 15 rare earths at pH 3, wash, acid desorption, ICP-MS;
  three replicates. Profile = log10 of desorbed amount, centred over the 15 elements.
- 353 of 600 orthologs pass quality control (total captured metal ≥ 20 ppb, replicate RSD ≤ 50%, two or more
  complete replicates). Two families hold 288 of them (Nitrobacteraceae 153, Methylobacteriaceae 135).
- My centred profile correlates 0.956 (median per protein) with the workbook's normalized logD; the two differ by
  a median 0.086 log10 units per value, so the workbook applies a normalization I did not reproduce.
- Average profile (log10): La +0.14, Ce +0.35, Pr +0.41, Nd +0.43, Sm +0.44, Eu +0.37, Gd +0.21, Tb +0.11,
  Dy −0.04, Ho −0.21, Er −0.34, Tm −0.43, Yb −0.48, Lu −0.56, Y −0.41.

## Step 0: this morning's frozen predictions

| Prediction (frozen blind) | Outcome |
|---|---|
| P1. Pooled gap between Gd and Tb below +15 mA | **Holds.** +10.0 mA [+10.0, +11.0]; per ortholog median +10, IQR +8 to +12 (small chelators: +31) |
| P2. Yttrium within −0.05 to +0.25 kcal/mol of a lanthanide of its radius; Dy or Ho closest | **Fails.** +0.378 kcal/mol [+0.370, +0.387]; closest are Tm (0.10) and Er (0.11) |

The profiles carry smooth higher-order curvature shared by all orthologs (a quartic leaves 0.031 log10), so the
single-break scan finds larger apparent breaks elsewhere (Nd|Sm −33 mA); none is a step at Gd|Tb. The revised,
non-blind expectation for yttrium (closest to Er or Tm) matches.

## Step 1: frozen first-look readings

| Reading | Outcome |
|---|---|
| A1. Reproducible variation (ICC over replicates) | **Pass.** Linear 0.73 [0.69, 0.77]; quadratic 0.51 [0.41, 0.59]; cubic 0.83 [0.80, 0.85]; yttrium offset 0.28 [0.18, 0.37] |
| A2. Profiles are smooth (cubic residual / replicate noise ≤ 2) | **Pass.** 1.47; residual RMS 0.051 log10. Differences between orthologs are smoother still: a quadratic leaves 0.94 of the noise level, a cubic 0.24 |
| A4. Sequence predicts profiles in unseen families | **Pass.** Nearest neighbour: skill +0.34 [+0.07, +0.55], cosine +0.56 [+0.35, +0.75] over 7 families. Language-model regression +0.13 [−0.26, +0.36]; loop-feature regression fails (−0.52) |

**A3. Size of natural differences** (log10 units; SD across orthologs after removing replicate noise):

| Contrast | Mean | SD | 5th–95th percentile |
|---|---|---|---|
| Lu − La | −0.71 | 0.364 | −1.23 to −0.02 |
| Dy − Nd | −0.46 | 0.122 | −0.66 to −0.24 |
| Y − Dy | −0.37 | 0.065 | −0.49 to −0.23 |
| Y − Ho | −0.20 | 0.038 | −0.27 to −0.11 |
| Tb − Gd | −0.10 | 0.034 | −0.16 to −0.05 |
| Dy − Tb | −0.14 | 0.031 | −0.20 to −0.07 |
| Nd − Pr | +0.01 | 0.016 | −0.01 to +0.04 |

- Nd over Dy ranges from a factor of 1.7 to 4.6 across orthologs. Seventeen of 353 prefer Lu to La; none
  prefers Dy to Nd.
- Genus explains 21% of the variance in Dy − Nd; most of it lies within genera.

**A5. How close must a neighbour be** (ten random folds; skill of the five-neighbour mean / language-model
regression): identity 0.95–1.0: +0.61 / +0.60 (n = 234); 0.85–0.95: +0.41 / +0.34 (77); 0.70–0.85: +0.40 / +0.46
(34); 0.50–0.70: −0.06 / −0.32 (8).

**A6. Sequence features** (frozen loop finder; within family): the linear tilt rises with acidic residues per
loop (+0.32 [+0.13, +0.41]) and with loops carrying Asp at position 9 (+0.39 [+0.08, +0.62]).

**A7. Capture against loop acidity:** Spearman +0.41 over all rows; see the clade check below.

## The three conditions stated in advance for "worth asking a laboratory"

1. **Natural lanmodulins differ reproducibly: met.**
2. **Within-branch contrasts or the yttrium offset vary by at least 0.1 log10: not met.** Nd − Pr 0.016,
   Dy − Tb 0.031, yttrium 0.04–0.07. Only cross-branch contrasts vary (Dy − Nd 0.122, Lu − La 0.364).
3. **Predictability at the identity of the environmental homologs: partly met.** Two have library relatives
   at 85–91% identity; the rest sit at 70% or below.

By the rule written beside those conditions, natural lanmodulin diversity is not worth mining for within-branch
or yttrium selectivity. The cross-branch tilt does vary, and the library already contains its extremes.

## Post hoc

- **Loop finder.** The frozen pattern required Asp at loop position 1 and agreed with the workbook's EF-hand
  count for 192 of 600 proteins; many lanmodulin loops start with Asn. A scored pattern allowing Asn agrees for
  403 of 600.
- **Tilt and loop acidity** (improved finder): Lu − La against acidic residues per loop +0.70 pooled, +0.56
  within family and within genus. The most Nd-over-Dy selective orthologs (factor 6–7) are Xanthobacter and
  Methylopila proteins with Asn at position 1 and four acidic residues per loop; the flattest are
  Methylobacterium proteins with five. The seal note of 9 November records that the article found only modest
  effects of swapping such loop residues, so this may not be causal.
- **Capture and loop acidity are confounded with clade.** The 92 three-EF-hand proteins (mostly Rhizobiaceae;
  3.35 acidic residues per loop) capture a median 4.7 ppb and 4% reach 20 ppb; the 508 four-EF-hand proteins
  capture 38.3 ppb and 79% reach 20 ppb. Within families and within genera the relation vanishes (−0.11
  [−0.17, +0.17]; −0.13 [−0.18, +0.02]). Hans-LanM, a picomolar lanmodulin with lower loop acidity, captures
  12 ppb here: capture at pH 3 penalizes less acidic sites and is not a test of binding at neutral pH.
- **Leave one genus out:** nearest neighbour +0.39 [+0.24, +0.50]; pooled over orthologs, five neighbours +0.58
  and language-model regression +0.51.

## Environmental homologs

The placement of unpublished environmental homologs against the library is kept in the private copy of this
report. In summary: two have library relatives at 85–91% identity and are predicted to have ordinary profiles; the
others lie at or beyond the edge of the predictable range.

## Deviations from PLAN.md

1. The frozen loop pattern performs poorly (above); loop-feature results as frozen are reported, and the
   improved finder is labelled post hoc.
2. Added after the frozen readings: smoothness of ortholog differences, leave-one-genus-out, the clade check on
   capture, the clade table, and the placement of environmental homologs.
3. Step 2 of the plan (locking a sequence model and scoring the 16 reserved proteins) was not run. No decision
   depends on it, and the 16 sealed proteins keep more value as an untouched test set.
4. `aligned identity × coverage` was used to place environmental sequences, because they include signal peptides.

## Limits

- One assay: competitive capture at pH 3 on a solid support. It measures selectivity among rare earths, with
  no calcium, and it under-reports proteins with less acidic loops.
- Two families make up 82% of the usable orthologs; leave-one-family-out has seven folds.
- The loop finder still disagrees with the workbook for a third of proteins.
- Function and metal binding of the environmental homologs remain UNVERIFIED.

## Files

`PLAN.md`, `PLAN.sha256`; `scripts/` (`extract_nonreserved.py`, `common.py`, `step0_frozen.py`,
`prep_sequences.py`, `first_look_a.py`, `first_look_b.py`, `capture_and_env.py`, `describe_clades.py`,
`make_figure.py`); `results/` (logs, `ortholog_components.tsv`, `ortholog_contrasts.tsv`, `loops_improved.tsv`,
summaries); `data/` (the 600-row extract, hashes in `EXTRACTION.json`);
`figures/lanmodulin_ortholog_first_look.png`.
