# Two-coefficient (radius-quadratic) selectivity model vs M1 — PLAN (frozen 2026-11-09, before any new fit)

**Question.** Does a two-coefficient model, which predicts the linear (a) and quadratic (b) radius coefficients from
ligand features and is trained on every measured lanthanide for the training ligands, predict held-out families better
than M1 (donor counts + charge ridge, family-held-out RMSE 0.854 kcal/mol, `diagnostics/ln_structure_increment_20261030`)?

Scope: login node, single thread, inline, seconds. No Slurm. Repository HEAD at freeze: 96952f236dfef4ffbce322a7d1eaf42388d13ba0
(dirty working tree from other sessions preserved; nothing committed, stashed or edited).
Code and this PLAN: `diagnostics/ln_two_coefficient_20261109/`. NIST-derived tables and logs (untracked):
`workspaces/ln_two_coefficient_20261109/`.

## Reserved / sealed material — never opened or used
- Reserved SpyCI-LAMBS (Spicy-Lams) assay outcomes and the 16 reserved-panel proteins:
  `diagnostics/spicy_lams_inventory_20260923/` and any `workspaces/spicy_lams*` path.
- Reserved PLM outcomes: `diagnostics/plm_*`, `diagnostics/nikasha_plm_*`, `workspaces/plm_*`, `workspaces/nikasha_plm*`.
- Treated as reserved out of caution (role unclear): `workspaces/ln_structure_increment_20261030/exclusions.tsv`
  (holds the TTHA/TMDTA and other exclusion rows). Not opened, not hashed. TTHA and TMDTA records are never used
  for fitting or scoring; NIST raw rows are filtered to the 115 admitted ligand ids before any value is parsed.
- `lit_untrusted/` directories: not read.
- Protein/peptide series (`protein_ln_affinities.tsv`, `lady_series_*`): not used.

## Inputs (only these)
- `workspaces/ln_structure_increment_20261030/primary.tsv` (115 primary La/Dy records, family, record_id/condition code,
  annotation code) and `.../features.tsv` (M1 feature columns). SHA256 of both checked against
  `diagnostics/ln_structure_increment_20261030/FREEZE.json` pins.
- `diagnostics/ln_structure_increment_20261030/{evaluate.py,curate.py,structure.py}` imported for `nested`, `ALPHAS`,
  `M1` (code only; their `main()` is not run).
- `workspaces/ln_selectivity_calibration_20261027/nist46_dat/` parsed with a copy of `curate.raw_records()` whose metal
  filter is widened from {La, Dy} to the 14 Ln(III) of the radius table (Pm excluded: no radius in the frozen table).
- Lig-id columns only of `workspaces/ln_selectivity_calibration_20261027/{features.tsv,pairs.tsv}` (R2-plus count).

## Data for the multi-metal fits (rule fixed now)
For each of the 115 ligands, take the condition code of its primary record (record_id suffix `_C<code>`). At that
condition, a metal M is admitted iff: exactly one native value; value is a plain number (same regex as curate.py);
annotation code equals the primary record's annotation code (the curate.py La/Dy rule). Target per admitted metal
y(L,M) = −KC·[log K1(M) − log K1(La)], KC = 1.3642470130342386 (curate.py). Admitted La/Dy values must equal
primary.tsv exactly (assertion). Other conditions of the same ligand are not used (they are linked secondary
records in the 30 Oct design, not independent examples). Counts of admitted and dropped metal values are reported.

## Positive control (Step 0; stop if it fails)
Rebuild M1 with the imported `nested()` on primary.tsv ⋈ features.tsv (evaluate.py merge), family = outer fold,
alpha grid [0.01, 0.1, 1, 10, 100], inner leave-one-family-out ligand-weighted MSE, StandardScaler inside folds.
Pass iff family-held-out RMSE is within 0.005 of 0.854.

## Model R2
y(L, M) = a(x_L)·φ1(M) + b(x_L)·φ2(M) + u_M, with Δr_M = r(M) − r(La) (Shannon CN9, table in
`diagnostics/ln_lattice_strain_20261109/fit_lattice_strain.py`), radius basis normalised to Dy so the single-metal case
maps exactly onto M1: φ1 = Δr/Δr_Dy, φ2 = (Δr/Δr_Dy)². (a, b are therefore Dy-scaled; a_Å = a/Δr_Dy, b_Å = b/Δr_Dy².)
a(x) = α0 + α·z, b(x) = β0 + β·z, z = M1 feature vector standardised with a scaler fitted on the training-fold ligands,
one row per ligand. u_M one unpenalised offset per non-La metal, u_La = 0. With u free, α0 and β0 are absorbed into u
(not separately identified) and are dropped; in the no-offset variant α0, β0 are kept, unpenalised.
Fit: least squares, one row per admitted (ligand, metal≠La) point, unweighted, plus ridge penalty λ·(‖α‖²+‖β‖²).
Grid λ ∈ 2×[0.01, 0.1, 1, 10, 100] (at Dy, φ1 = φ2 = 1, so penalty λ on (α, β) equals sklearn alpha = λ/2 on their
sum; this makes R2-dy algebraically identical to M1). Nested selection: inner leave-one-family-out within the outer
training families; criterion = ligand-weighted MSE of Dy−La on the inner held-out ligands (same criterion as M1);
lowest λ wins ties. Outer folds = the 11 families of the 30 Oct splits.

Variants, all scored on the same held-out folds:
- R2-dy: training rows = the 115 ligands' Dy−La only. Check: |RMSE − 0.854|.
- R2-all: every admitted Ln at the primary condition for the training-fold ligands.
- R2-plus: R2-all plus ligands that have feature vectors (lig_id in either features table) but lack La or Dy.
  Count reported; if zero, skipped.
- R2-all-nou: R2-all with u_M ≡ 0 (α0, β0 unpenalised).

## Scores
- Primary: family-held-out RMSE of Dy−La on the 115 ligands, prediction a·φ1(Dy) + b·φ2(Dy) + u_Dy.
- Secondary: held-out ligands with complete series (all 13 non-La metals admitted at the primary condition): RMSE over
  their 13 metals. Comparator "M1 → SVD shape 1": V1 = first right singular vector of the uncentred matrix of the
  training-fold complete-series ligands (13 metals); series = ŷ_Dy(M1)·V1/V1[Dy].
- Uncertainty: paired bootstrap of RMSE(M1) − RMSE(variant), 2,000 resamples, numpy default_rng(42), 95% percentile
  intervals, (i) over ligands and (ii) over families (whole families resampled). Same for the secondary score.

## Frozen readings (primary score, R2-all)
- "Strong": RMSE ≤ 0.60.
- "Supported": RMSE(M1) − RMSE(R2-all) ≥ 0.05 with the ligand-level interval excluding 0.
- "Not supported": anything else.

## Descriptive (no reading)
- Per-family RMSE, M1 and R2-all.
- Residuals (obs − pred) under M1 and R2-all for TTHA, CDTA, TMDTA, DTPA-bis(methylamide),
  oxolane-2,3,4,5-tetracarboxylate, DTPA, EDTA where in the 115 (TTHA and TMDTA are excluded upstream; reported absent).
- Per-ligand F2 fit y = a·Δr + b·Δr² (Å units, through origin, admitted metals) for ligands with ≥4 metals including La;
  family-held-out prediction of a and of b from z with the same nested ridge (sklearn, M1 grid);
  R² = 1 − Σ(obs−pred)²/Σ(obs−mean)².
- Selected λ per fold; u_M estimates.
