# Lattice-strain test of the La–Lu two-shape basis — PLAN (frozen 2026-11-09, before any fit)

**Question.** Is the empirical two-shape SVD basis of NIST SRD 46 ΔΔG(M−La) the linear and quadratic terms of a
lattice-strain parabola in ionic radius (Brice 1975; Blundy & Wood 1994)? Does a physically constrained
two-parameter form in radius describe the series as well as the SVD basis does?

Scope: descriptive fitting of existing tables on the biotite login node, single thread, inline. No Slurm, no new
chemistry, no reserved/sealed/holdout tables.

## Inputs (only these)
- `workspaces/ln_surrogate_20261027/surrogate_nist.tsv` parsed exactly as `diagnostics/lady_series_20261028/series_rank.py`
  (47 ligands with complete La..Lu series; ΔΔG(M−La) = −1.3643·[log K1(M) − log K1(La)] kcal/mol, 13 metals).
- `workspaces/ln_selectivity_calibration_20261027/protein_ln_affinities.tsv`, datasets built exactly as
  `diagnostics/lady_series_20261028/protein_series_shapes.py` (code copied, script itself not executed, so its output
  table is left untouched).
- `workspaces/lady_series_20261028/protein_series_shapes.tsv` (existing two-shape residuals, read only).
- `workspaces/ln_selectivity_calibration_20261027/features.tsv` (donor counts, charge, carboxylates, denticity).

## Step 0 — positive controls (stop if either fails)
- PC1: rerun `series_rank.py` unchanged: SVD variance 0.981 / 0.013 / 0.003; LOO RMS 0.28 (Dy) and 0.12 (Nd+Dy+Lu).
- PC2: my copy of the protein dataset builder reproduces `res2` in `protein_series_shapes.tsv` within 0.002 for all 18 rows.

## Radii (Shannon 1976 effective radii, Å, Ln³⁺)
CN9: La 1.216, Ce 1.196, Pr 1.179, Nd 1.163, Sm 1.132, Eu 1.120, Gd 1.107, Tb 1.095, Dy 1.083, Ho 1.072, Er 1.062, Tm 1.052, Yb 1.042, Lu 1.032.
CN8: La 1.160, Ce 1.143, Pr 1.126, Nd 1.109, Sm 1.079, Eu 1.066, Gd 1.053, Tb 1.040, Dy 1.027, Ho 1.015, Er 1.004, Tm 0.994, Yb 0.985, Lu 0.977.
Verification (done before fitting): all 28 values equal the Shannon 1976 entries in the `mendeleev` package database
(lnsel env, `element(X).ionic_radii`, charge 3, CN VIII / IX). CN9 contractions La→Nd 0.053, La→Dy 0.133, La→Lu 0.184
match the lady_series README.

## Forms (per ligand, y(M) = ΔΔG(M−La), Δr = r(M) − r(La); every NIST form passes through y(La) = 0)
- **F1** SVD two-shape basis: shapes V1, V2 from the uncentred SVD of the other 46 ligands (leave-one-ligand-out);
  y = c1·V1 + c2·V2. 2 parameters.
- **F2** y = a·Δr + b·Δr², CN9 radii. 2 parameters.
- **F3** same with CN8 radii. 2 parameters.
- **F4** Brice: y = S(r_M) − S(r_La), S(r) = 1.809·E·[(r0/2)(r−r0)² + (1/3)(r−r0)³], E in GPa, r in Å, y in kcal/mol, CN9 radii.
  Unit factor: 4π·N_A·(1e9 Pa)·(1e-30 m³) = 7567.6 J/mol = 1.8087 kcal/mol per GPa·Å³ (checked).
  Fit: for each r0, E is linear (closed-form least squares, sign free); r0 profiled on a grid 0.300–3.000 Å step 0.001,
  then refined with a bounded scalar minimiser within ±0.002 Å of the best grid point. 2 parameters.
  Fits with r0 within 0.005 Å of a grid bound are flagged "boundary (degenerate)". Unphysical: r0 outside 0.8–1.4 Å or E ≤ 0.
  Secondary variant **F4c** (reported, no reading): E ≥ 0 and r0 constrained to 0.8–1.4 Å.
- **F5** y = a·z + b·z², z = Z − 57. 2 parameters.
- **F6** continuous piecewise-linear in Δr (CN9), break at Gd: y = a·Δr + b·min(0, Δr − Δr_Gd). 2 parameters.
- **F7** y = a·Δr (CN9). 1 parameter.
- **F8** y = c1·V1 (LOO SVD shape 1). 1 parameter.
- Reference (existing method, reported beside F1/F8): the `series_rank.py` LOO regression map from the fit metals to
  all 13 metals (F1r for Nd+Dy+Lu, F8r for Dy).

## Scores (kcal/mol)
- **S1** in-sample RMS residual over all 47 × 13 points (F1/F8 use LOO shapes; radius forms fitted per ligand on all 13 metals).
- **S2** fit each 2-parameter form on Nd, Dy, Lu only; predict the other 10 metals; RMS over 47 × 10 held-out points.
- **S3** fit each 1-parameter form on Dy only; predict the other 12 metals; RMS over 47 × 12 held-out points.
- Secondary convention, reported for comparability with `series_rank.py`: S2/S3 RMS over all 13 metals (fit metals included).
- **Uncertainty**: bootstrap over ligands (resample the per-ligand squared-error vectors), 2,000 resamples,
  `numpy.random.default_rng(42)`, 95% percentile intervals. Differences are paired (same resampled ligands).
  F1/F8 LOO shapes are computed once on the full set and their LOO errors resampled.

## Frozen readings
- **Adequate**: best of F2/F3/F4 has S2 ≤ 0.20 kcal/mol.
- **Equivalent to SVD**: that S2 is within +0.03 of F1's S2 (primary: F1 as defined above; also reported against F1r).
- **Smooth bend**: F2 beats F6 on S1 with the bootstrap interval of S1(F6) − S1(F2) excluding 0 (positive).
  **Kink**: the reverse (interval of S1(F6) − S1(F2) entirely negative). Otherwise **undetermined**.

## Descriptive results (no gate)
1. F2 (CN9) in-sample residuals: mean per metal across ligands with SE = SD/√47; |t| per metal; share of residual
   sum of squares carried by the per-metal means; sign-flip test (flip each ligand's residual vector, 2,000, seed 42)
   on Σ_M mean(M)². Pattern described against a Gd break and the tetrad groups (La–Nd, Pm–Gd, Gd–Ho, Er–Lu).
   Universal per-metal offset c(M) = mean of the other 46 ligands' F2 residuals (LOO): S1 and S2 recomputed with
   y − c(M) fitted and c(M) added back.
2. F4: distributions (median, IQR) of r0 and E; counts of boundary, unphysical (r0 ∉ [0.8, 1.4] or E ≤ 0).
   "Monotonic ligand" := F2 vertex −a/(2b) outside the La–Lu Δr range (or b = 0). If F4 is degenerate for monotonic
   ligands, (a, b) from F2 serve as the stable parameters.
3. Spearman ρ (bootstrap 95% intervals, 2,000, seed 42) of F2 a and b (and physical F4 r0, E) with: total donor
   groups (sum of the 12 donor-group columns carboxyl … hydroxamate), ligand charge `q_lig`, `carboxyl`, denticity
   `nlig_La_cn8` (the column `fit_surrogate.py` uses as denticity), over ligands present in features.tsv.
4. Protein/peptide: the 18 datasets of `protein_series_shapes.py`, same metal subsets. Each fitted as
   y = c + f(r) (c free, as the existing script does), for F2–F7 (F4 as c + S(r)); residual RMS beside the existing
   res1/res2; shortcut (fit on La, Nd, Dy, Lu where all present and n > 4, predict the rest) beside `short`.
   NIST reference for F2 on the same subsets: median per-ligand residual of the same fit.

## Outputs
- Code, PLAN, REPORT (≤150 lines), DEVIATIONS (only if needed): `diagnostics/ln_lattice_strain_20261109/`.
- Data tables (NIST-derived, untracked): `workspaces/ln_lattice_strain_20261109/`.
