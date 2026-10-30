# NIST SRD 46 lanthanide selectivity data and the donor-count baseline (2026-10-27; audit 2026-10-29)

Empirical La/Dy selectivity from the Martell–Smith critically selected constants (NIST SRD 46 v8). It is the baseline every computed method has been compared against, and the data source for the predictor proposed in `../../LANTHANIDE_NIST_PREDICTOR_REVIEW_2026-10-29.md`.

Convention: dlogK = log K1(Dy) − log K1(La); ΔΔG(Dy−La) = −1.364·dlogK kcal/mol at 25 °C, so negative ΔΔG means Dy is favoured.

## Pipeline

The scripts run inside `workspaces/ln_selectivity_calibration_20261027/` (not tracked) with `/groups/banfield/users/jwestrob/envs/lnsel/bin/python`. All steps are single-threaded and take seconds, apart from the PubChem queries.

| Step | Output | Content |
|---|---|---|
| `extract_ln.py nist46_dat` | `ln_ML.tsv` | Every Ln(III) 1:1 log K, ΔH and ΔS row: 6,033 rows; 338 ligands have a 25 °C log K |
| `build_pairs.py` | `pairs.tsv` | 171 ligands with La and Dy log K1 at 25 °C and the same ionic strength (the one closest to 0.1 M) |
| `fetch_smiles.py` | `smiles.tsv` | PubChem SMILES by name, kept only if the formula matches NIST's: 113 of 171 |
| `fit_donors.py` | `features.tsv`, `fit_report.txt` | 12 RDKit donor-group counts; ridge regression scored by leave-one-out |
| `data_audit_20261029.py` | `DATA_AUDIT_20261029.txt` | The numbers behind the predictor proposal (runs from any directory) |

The NIST data files come from the free SRD 46 installer (https://www.nist.gov/srd/nist46), decoded without running it by `../lady_relaxed_exchange_20261024/tools/nist46.py`.

**The NIST-derived tables are not committed.** This repository is public, so the scripts regenerate the tables from the installer rather than redistributing the compilation. Documents quote individual values.

## Result (2026-10-27, `fit_report.txt`)

The set is n = 101 ligands (I ≤ 0.5, formula-checked SMILES), with sd(ΔΔG) 1.12 kcal/mol. Leave-one-out RMSE:

| Model | RMSE (kcal/mol) |
|---|---|
| mean | 1.14 |
| log K_La alone (stronger binders are more Dy-selective) | 0.81 |
| 12 donor-group counts | 0.85 |
| donor-group counts + log K_La | 0.81 |

log K_La is measured, not predicted, so it is not available for a new ligand.

## Audit (2026-10-29, `DATA_AUDIT_20261029.txt`)

- **The error is concentrated.** The 7 largest residuals carry 63% of the donor model's squared error, and the other 94 ligands are fit to 0.54 kcal/mol.
  - The seven are TTHA (−4.2: predicted far too Dy-selective), CDTA (+3.2), TMDTA (+2.35), DTPA-BMA (+2.2), oxolane-2,3,4,5-tetracarboxylate (+2.0), DTPA (+1.5) and EDTA (+1.4).
  - All are polycarboxylates whose selectivity is set by the backbone and the number of donors, which the counts do not encode.
- **Same donors, different backbone.** The 15 diamine-N,N,N′,N′-tetraacetates span dlogK 0.99–4.16, that is 4.3 kcal/mol.
  - The donor model predicts 1.8–2.0 for the four of them that it contains.
  - The other 11 are missing: 9 because the PubChem lookup failed, 2 because they were measured at I = 1.0.
- **Data the current model leaves out:**
  - 58 La/Dy pairs whose PubChem lookup failed. These are mostly backbone variants: cycloalkylene isomers, longer linkers, propionate arms, DTPA bis-amides, aza-crown acetates.
  - 12 pairs measured at I > 0.5.
  - 37 ligands that lack La or Dy but have a light (La–Nd) and a heavy (Tb–Lu) Ln under the same conditions, so they can enter through the one-shape series model (`../lady_series_20261028/README.md`).
  - That gives 208 ligands in total.
- **Noise floor.** 33 ligands are listed at two or more ionic strengths. Their dlogK moves by a median of 0.09 log units, and the per-value sd is 0.29 kcal/mol, or 0.57 for strong chelators (log K_La > 10). This mixes real ionic-strength dependence with between-study scatter, so read it as the floor for pooled data.
