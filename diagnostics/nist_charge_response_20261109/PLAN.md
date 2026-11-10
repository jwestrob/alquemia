# NIST charge response: La³⁺ versus Ca²⁺ and other metals — PLAN (frozen 2026-11-09, before any fit)

**Question.** Does log K1(La) − log K1(Ca) (two ions of nearly equal size, Shannon CN8 1.16 and 1.12 Å, differing by
one unit of charge) form an experimental selectivity scale that donor-group counts predict for held-out ligand
families? Secondary: the same question for absolute log K(La), log K(Ca) and La-versus-other-metal differences.

Scope: biotite login node, single thread, inline, seconds. No Slurm, no GPU. Repository HEAD at freeze:
7b028c3b9cb59b3033bc5150d73ee6826c48b06d (dirty working tree from other sessions preserved; nothing committed,
stashed or edited). Code and this PLAN: `diagnostics/nist_charge_response_20261109/`. NIST-derived tables and logs
(untracked): `workspaces/nist_charge_response_20261109/`.

## Reserved / sealed material — never opened or used
- Reserved SpyCI-LAMBS assay outcomes and reserved-panel proteins: `diagnostics/spicy_lams_inventory_20260923/`,
  any `workspaces/spicy_lams*` path.
- Reserved PLM outcomes: `diagnostics/plm_*`, `diagnostics/nikasha_plm_*`, `workspaces/plm_*`, `workspaces/nikasha_plm*`.
- Treated as reserved out of caution: `workspaces/ln_structure_increment_20261030/exclusions.tsv` (not opened, not
  hashed; the FREEZE pin check covers the other 11 of 12 files).
- Every `lit_untrusted/` directory: not read.
- Protein/peptide tables (`protein_ln_affinities.tsv`, `lady_series_*`): not used.

## Extraction (already run, counts only inspected; `extract_pairs.py`, log `extract.log`)
Equilibrium ML = [ML]/[M][L], quantity K, T = 25 °C, from `nist46_dat/F8558_VALUES.DAT` and code tables.
For each ligand and partner metal X: if La and X share condition codes, the 30 Oct curate rule applies (exactly one
value per metal, plain number, equal annotation codes); among admissible shared conditions the one nearest I = 0.1 M
is primary (ties: lower I, then lower code). If shared conditions exist and none is admissible, the ligand is
excluded for X. If La and X share no condition, the pair with the smallest |I_La − I_X| among single plain values is
recorded (`nearest_I`, ΔI and annotation equality saved). Positive controls passed before this PLAN: the historical
build_pairs rule gives 171 La/Dy ligands (pairs.tsv 171); the strict rule reproduces all 115 primary La/Dy records
of the 30 Oct design exactly (record id and both values).

Feature vectors and families: the 115 ligands of `workspaces/ln_structure_increment_20261030/{primary,features}.tsv`
only (M1 columns, 30 Oct family labels, normalized identities). No new SMILES (no web lookups); no new families.

## Targets (log units; positive = La binds more strongly)
- **y1** = log K(La) − log K(Ca) (primary). y2 = log K(La); y3 = log K(Ca) (same pair records).
- Secondary: y4 = La − Sr; y5 = La − Al; y6 = La − Fe³⁺. Added before fitting from the database inventory:
  y7 = La − Mg, y8 = La − Ba (divalent size series), y9 = La − Lu (same charge, size only),
  y10 = La − Cd (Ca-sized soft divalent), y11 = La − Y.
- **Scored set (primary)** per target: featured ligands with a same-I pair, plus nearest-I pairs with |ΔI| ≤ 0.1 M.
  Rationale fixed now: Debye–Hückel predicts the La−Ca difference shifts by about 2|z_L|·A·f(I), which is
  ~0.2|z_L| log units between I = 0 and 0.1 M, so larger gaps are excluded from the primary set.
  Sensitivity sets: same-I only; all pairs (any ΔI). Both reported, no reading.
- A target is scored only if its primary set has n ≥ 25 ligands in ≥ 4 families; otherwise counts only.

## Models (all with the imported 30 Oct `nested()`: outer leave-one-family-out, inner leave-one-family-out,
ligand-weighted MSE, StandardScaler inside folds, alpha grid [0.01, 0.1, 1, 10, 100], lowest alpha wins ties)
- M0: training-fold mean. M1: 12 donor SMARTS counts + macrocycle + ligand charge (30 Oct `structure.M1`).
- D_anO: one feature, anionic O-donor groups = carboxyl + phosphonate + phenol + hydroxamate (sulfonate is a
  spectator, as in the 30 Oct geometry rule).
- D_tot: one feature, total donor groups = carboxyl + amine_N + ether_O + alcohol + phenol + phosphonate +
  pyridine_N + amide_O + carbonyl_O + hydroxamate (aminoacid_N omitted: it re-counts amine N; sulfonate omitted).

## Positive control (Step 0; stop if it fails)
M1 on the 115 La/Dy primary records (kcal/mol) must give family-held-out RMSE within 0.005 of 0.854.

## Scores
Family-held-out RMSE and R² = 1 − SSE/Σ(y − ȳ)² (ȳ = mean of the evaluated set) per model. Paired bootstrap,
2,000 resamples, 95% percentile intervals: (i) over ligands, numpy default_rng(42) integers; (ii) over families,
default_rng(42) choice of whole families. Intervals for each model's RMSE and R², and for RMSE(M0) − RMSE(M),
M ∈ {M1, D_anO, D_tot}, and RMSE(D_anO) − RMSE(M1).

Noise floor: ligands with admissible same-condition La/Ca pairs at ≥ 2 ionic strengths; pooled within-ligand SD of
y1 (degrees of freedom n_i − 1) and median within-ligand range. It mixes real ionic-strength dependence with
between-study scatter.

## Frozen reading (y1, primary set, model M1)
- "strong scale": held-out R² ≥ 0.6.
- "usable": R² ≥ 0.3 and the 95% interval on RMSE(M0) − RMSE(M1) excludes 0 under BOTH the ligand and the family
  bootstrap (stricter choice made here because the lead's wording names one interval without specifying which).
- "not supported": anything else.

## Descriptive (no reading)
- M1 ridge coefficients for y1 in raw units (log units per group = coef/scale), alpha chosen by leave-one-family-out
  on the whole primary set with the same criterion; 2,000 ligand-bootstrap refits at that alpha (default_rng(42)),
  95% percentile intervals; family-bootstrap intervals alongside. Groups are partially collinear (aminoacid_N ⊂ amine_N).
- y1 by carboxylate count: n, median, min, max.
- Ten largest and ten smallest y1 with nonzero donor groups and charge.
- Spearman ρ between y1 and the lattice-strain linear coefficient a2
  (`workspaces/ln_lattice_strain_20261109/nist_per_ligand_fits.tsv`) on shared ligands, ligand-bootstrap interval.
- Penalty-grid edge: selected alpha per outer fold. If any fold selects 0.01 or 100 for M1, rerun with
  [0.001 … 1000] and report as OUTSIDE the frozen design.
