# Addendum: confirmation on ligands not yet used (frozen 2026-11-10, before those series were read)

Written by the lead agent after the frozen readings of PLAN.md were scored and after the post hoc analyses on
the 47 complete-series ligands (`scripts/run_followup.py`, `fit_cubic_coordinate.py`, `fit_step.py`,
`explore_groups.py`, `yttrium_branch.py`). Those analyses suggested three things that PLAN.md did not anticipate:

1. Each ligand needs a cubic in radius (four numbers).
2. With a cubic per ligand, the shared structure is a single extra gap between Gd and Tb: the light branch
   (La–Gd) sits +31.0 mA from the heavy branch (Tb–Lu) on the Shannon CN9 scale (bootstrap 95% +18.5 to +61.5;
   best break position in 945 of 1,000 resamples). Cross-validated error on the 47: cubic 0.1584, cubic with the
   gap 0.1382, cubic in a free 12-number coordinate 0.1356 kcal/mol.
3. Yttrium (Shannon CN9 radius 1.075 A) sits on the heavy branch for oxygen-only ligands (fraction of the gap
   0.18, 95% −0.10 to 0.52) and on the light side for nitrogen-containing ligands (1.14, 0.68 to 2.50).

All three came after seeing results, so they are tested here on data that no fit has touched: the 24 ligands of
`data/surrogate_nist.tsv` with 10 to 13 of the 14 lanthanides (only their metal counts have been read). 17 of
them list Y.

## Fixed inputs (no refitting on the 24)

- Gap: +0.0310 A added to the Shannon CN9 radius of La…Gd (model KS).
- Free coordinate: `results/coordinate_cubic.tsv`, column `rho_KC`
  (sha256 401876b2b3aa76ea0263edb1b57fead5a4f8db792131517b7a8f5e26f63b8bff) (model KC).
- Shared offset: mean series of the 47 ligands (model KO).
- N atoms counted from the molecular formula.

## Tests

For each of the 24 ligands, leave-one-metal-out over its available lanthanides: fit the ligand's own four numbers
on the rest, predict the held-out metal. RMS per ligand; bootstrap over the 24 ligands (2,000 resamples, seed 42),
paired differences.

- **C1 (primary). The gap.** Confirmed if RMS(K) − RMS(KS) is positive with its 95% interval excluding 0.
- **C2. The free coordinate.** RMS(K) − RMS(KC), same criterion; RMS(KS) − RMS(KC) reported.
- **C3. Independent estimate.** The gap fitted on the 24 ligands alone, with a bootstrap interval, and the best
  single break position among the twelve neighbour pairs. Replicates if the interval excludes 0 and Gd|Tb is the
  best position.
- **C4. Offset.** RMS(K) − RMS(KO); the 47-ligand result predicts no gain.
- **C5. Yttrium.** For the ligands that list Y, predict y(Y) from the ligand's lanthanides alone (cubic in the
  two-branch coordinate, all available lanthanides) at radius 1.075 on the heavy branch and at 1.075 + 0.031 on
  the light branch. Replicates if RMS(light) < RMS(heavy) for N-containing ligands and RMS(heavy) < RMS(light)
  for O-only ligands; paired intervals reported. Also the best single fraction per class with bootstrap interval.
- **C6. Proteins and peptides.** The 12 datasets with at least 9 lanthanides: leave-one-metal-out RMS under K and
  KS (gap fixed at +31.0 mA). Transfers if KS is better in a two-sided sign test at p < 0.05.

Quadratic counterparts (Q against the quadratic with the gap) are reported beside each cubic result.
