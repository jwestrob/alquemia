# Lattice-strain test of the La–Lu series shapes — 9 November 2026

**Result.** A two-parameter quadratic in ionic radius describes each ligand's La–Lu selectivity series
about as well as the two-shape SVD basis. Both frozen readings pass, and the bend is smooth.

All numbers are in `workspaces/ln_lattice_strain_20261109/{run.log,summary.json}`. This report was
written by the lead agent from those files (see DEVIATIONS.md).

## Frozen readings

`PLAN.md` (sha256 `4aee7648…f106`) was written before any fit and is unchanged.

| Reading | Criterion | Outcome |
|---|---|---|
| Adequate | best radius form S2 ≤ 0.20 kcal/mol | **Pass.** Brice form (F4) S2 = 0.183 [0.151, 0.214] |
| Equivalent to SVD | S2 within +0.03 of the two-parameter SVD fit (F1) | **Pass.** F4 − F1 = +0.019 [−0.002, +0.039] |
| Smooth bend or kink at Gd | S1(F6) − S1(F2), interval excluding 0 | **Smooth.** +0.032 [+0.015, +0.050] |

Positive control: `series_rank.py` reproduces 0.981 / 0.013 / 0.003, Dy 0.28 and Nd+Dy+Lu 0.12.

## Scores

kcal/mol, 47 ligands with complete series. S1 = in-sample. S2 = fitted on Nd, Dy, Lu and scored on the
other metals. S3 = fitted on Dy and scored on the others. Intervals: bootstrap over ligands, 95%.

| Form | Parameters | S1 | S2 | S3 |
|---|---|---|---|---|
| F1 SVD two-shape | 2 | 0.138 | 0.164 [0.132, 0.195] | |
| F2 radius quadratic, CN9 | 2 | 0.164 | 0.186 [0.153, 0.217] | |
| F3 radius quadratic, CN8 | 2 | 0.173 | 0.195 [0.161, 0.229] | |
| F4 Brice (r0, E) | 2 | 0.161 | 0.183 [0.151, 0.214] | |
| F5 quadratic in Z | 2 | 0.200 | 0.218 [0.180, 0.256] | |
| F6 kink at Gd | 2 | 0.196 | 0.224 [0.182, 0.264] | |
| F7 linear in radius | 1 | 0.311 | | 0.352 |
| F8 SVD shape 1 | 1 | 0.242 | | 0.291 |

The existing regression map (39 coefficients shared across ligands) scores 0.140 on the held-out
metals and 0.123 over all 13, matching the README's 0.12. The radius forms sit +0.043 above it on the
held-out metals; it is a different kind of model from a two-parameter per-ligand fit.

CN9 radii fit better than CN8 (S2 difference −0.010 [−0.014, −0.005]), and radius fits better than
atomic number (−0.032 [−0.047, −0.018]).

## What the two shapes are

- SVD shape 1 lies 99.7% inside the plane spanned by Δr and Δr². Shape 2 lies 95.7% inside it.
- After the quadratic, every ligand shares one residual pattern by metal, about 0.1 kcal/mol in size
  (sign-flip p = 0.0005; 28% of the residual variance):

| Ce | Pr | Nd | Sm | Eu | Gd | Tb | Dy | Ho | Er | Tm | Yb | Lu |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| −0.10 | −0.12 | −0.10 | −0.07 | −0.03 | +0.09 | +0.07 | +0.08 | +0.09 | +0.08 | −0.01 | −0.09 | −0.10 |

  Standard errors are 0.01–0.03. A cubic term explains R² 0.90 of this pattern, a step at Gd 0.60, and
  the two together 0.95.
- Quadratic plus that universal per-metal offset (learned leave-one-ligand-out) gives S2 = 0.165,
  equal to the two-parameter SVD fit (0.164).

## Parameters

- **Linear coefficient a** (kcal/mol/Å) follows ligand charge and carboxylate count: ρ = 0.78
  [0.60, 0.89]. Donor total: 0.70 [0.50, 0.84].
- **Quadratic coefficient b** follows denticity weakly: ρ = 0.36 [0.02, 0.65].
- **Brice (r0, E)** is physical for all 22 ligands whose optimum lies inside the series and for 8 of
  the 25 monotonic ligands. Over those 30: r0 median 1.061 Å (IQR 1.022–1.106), E median 35.9 GPa
  (IQR 32.4–51.5, range 12.5–215). E follows denticity (ρ = 0.63 [0.28, 0.86]); r0 falls with charge
  (ρ = −0.44 [−0.73, −0.05]).
- For monotonic series r0 is unconstrained (profile width 1.7 Å against 0.04 Å for the others), so
  (a, b) are the stable parameters across all ligands.
- Worked example, acetate: a = 7.93, b = 46.5. Dy−La = −0.133·a + 0.0177·b = −1.05 + 0.82 =
  −0.23 kcal/mol fitted (NIST −0.07). A linear gain is offset by a quadratic penalty of similar size.

## Protein and peptide series (18 datasets)

- The radius quadratic fits 16 of 18 at or below the stored two-shape residual (median 0.128 against
  0.156). The stored residuals reproduce to 0.0005.
- Shortcut from La, Nd, Dy and Lu (12 datasets): median 0.158 against 0.185.
- Poor fits remain for GGR (0.36) and one Deblonde 2020 HisL series (0.81, from 1.31 under two shapes).

## Interpretation (inference)

The two-shape basis is a radius quadratic plus a ligand-independent per-metal pattern. Each site
reduces to a linear radius response, which follows charge, and a curvature, which follows denticity.
The origin of the per-metal pattern is unverified; the aqua-ion reference, shared by every ligand, and
the radius scale itself are the two candidates.

## Limits

- Radii were checked against one database (`mendeleev`, Shannon 1976); all 28 match.
- Charge and denticity came from `surrogate_nist.tsv` (MACE-built complexes), because `features.tsv`
  carries neither.
- Correlations with features are in-sample and descriptive (n = 47, or 30 for Brice parameters).

## Files

- `diagnostics/ln_lattice_strain_20261109/`: `PLAN.md`, `fit_lattice_strain.py`, this report, `DEVIATIONS.md`.
- `workspaces/ln_lattice_strain_20261109/`: `run.log`, `summary.json`, `nist_per_ligand_fits.tsv`,
  `protein_lattice_strain.tsv` (untracked; NIST-derived).
