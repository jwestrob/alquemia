# Two-coefficient predictor test — 9 November 2026

**Result: not supported.** A model that predicts the linear and quadratic radius coefficients from
ligand features, trained on every measured lanthanide, predicts held-out Dy−La worse than the
donor-count baseline M1: RMSE 0.914 against 0.854 kcal/mol.

Authorized by Jacob on 9 November; this reopens predictor work closed on 30 October for one bounded
test. Written by the lead agent from the analysis subagent's reply and the saved tables (see
DEVIATIONS.md). `PLAN.md` sha256 `6aca312d…0693` was written before any new fit.

## Positive control

M1 rebuilt with the 30 October `nested()` code on `primary.tsv` joined to `features.tsv` gives 0.8545
(target 0.854). FREEZE pins match for 11 of 12 files; the twelfth, `exclusions.tsv`, was left unopened
as possibly reserved.

## Data

Every lanthanide measured at each ligand's primary condition, under the `curate.py` rules: 1,302
non-La points over all 115 ligands, median 13 metals per ligand, 57 ligands with complete series.
Dropped: 37 censored or flagged values, 9 with an annotation differing from the primary record, 3
duplicates, and 16 Pm rows (no radius).

## Primary score

Dy−La RMSE, kcal/mol, n = 115, 11 families, family-held-out. Δ = M1 − variant (positive = variant
better). Paired bootstrap, 95%.

| Variant | RMSE | Δ | Ligand interval | Family interval |
|---|---|---|---|---|
| R2-dy (Dy−La only) | 0.854 | 0.000 | [0, 0] | [0, 0] |
| **R2-all (all metals)** | **0.914** | **−0.060** | [−0.177, +0.055] | [−0.176, +0.056] |
| R2-all-nou (no per-metal offsets) | 0.939 | −0.085 | [−0.208, +0.043] | [−0.204, +0.045] |
| R2-plus | skipped | | | |

- R2-dy matches M1's predictions to 4e−14, so the implementation is exact.
- R2-plus was skipped: all 122 ligands with feature vectors already have both La and Dy.
- Frozen reading: "supported" needed Δ ≥ 0.05 with the ligand interval excluding 0. **Not supported.**

## Secondary score

57 complete-series held-out ligands, 8 families, RMSE over 13 metals. The route "M1 then SVD shape 1"
scores 0.608.

| Variant | RMSE | Δ | Ligand interval |
|---|---|---|---|
| R2-all | 0.724 | −0.116 | [−0.286, +0.061] |
| R2-all-nou | 0.727 | −0.120 | [−0.289, +0.057] |

## Descriptive results

- **Coefficients from features** (family-held-out, n = 115): a reaches R² +0.41 (Spearman +0.53); b
  reaches R² −0.05. The curvature is unpredictable from these features.
- **Per family, M1 → R2-all.** Better: aliphatic carboxylates 0.621 → 0.525, aromatic carboxylates
  0.980 → 0.885, acyclic carbonyl 1.364 → 0.749 (n = 1). Worse: hydroxy carboxylates 0.786 → 1.055,
  ether carboxylates 1.006 → 1.313, aminocarboxylates 1.056 → 1.121, macrocycles 0.589 → 0.834,
  nitroso phenols 0.655 → 0.794.
- **Named ligands, residual M1 → R2-all.** CDTA −3.20 → −2.11; DTPA-BMA −2.47 → −1.83; DTPA −1.51 →
  −0.20; EDTA −1.53 → −0.45; oxolane tetracarboxylate −1.64 → −2.14. TTHA and TMDTA are outside the
  115.
- **Penalty.** Every outer fold of every model, M1 included, selected the top of the grid (λ = 200).
  The frozen grid edge binds.

## Interpretation (inference)

The extra metals pin down each training ligand's curvature, but the features cannot predict curvature
for a new ligand, and the shared shrinkage gives Dy−La little benefit. The aminopolycarboxylates with
large residuals improve; losses in the hydroxy and ether carboxylates outweigh that. With the offsets
free, the fitted u_M absorb the mean curve (about −1.2 at Dy); they are a different quantity from the
0.1 kcal/mol universal residual in the lattice-strain report.

## Limits

- Whether `exclusions.tsv` is reserved was left unresolved; it was treated as reserved.
- A wider penalty grid is untested. The edge binds for M1 too.

## Reserved material left unopened

`diagnostics/spicy_lams_inventory_20260923/` and `workspaces/spicy_lams*`; `diagnostics/plm_*`,
`diagnostics/nikasha_plm_*`, `workspaces/plm_*`; `workspaces/ln_structure_increment_20261030/exclusions.tsv`;
all `lit_untrusted/` directories; the protein series.

## Files

- `diagnostics/ln_two_coefficient_20261109/`: `PLAN.md`, `fit_two_coefficient.py`.
- `workspaces/ln_two_coefficient_20261109/` (untracked): `run.log`, `summary.json`, `predictions.tsv`,
  `per_family.tsv`, `per_ligand_ab.tsv`, `multi_metal_points.tsv`.
