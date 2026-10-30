# Lanthanide-series extension of the La/Dy panel (2026-10-28)

## Why three metals are enough

Across the 47 NIST 46 ligands with complete La–Lu series, selectivity relative to La is close to one-dimensional. An uncentred SVD of ΔΔG(M−La) over 13 metals:
- Shape 1, a steady rise with contraction: 98.1% of variance.
- Shape 2, a bend at Gd/Tb: 1.3%.
- Shape 3: 0.3%.

The table shows how well a few computed metals reconstruct all 13 (leave-one-ligand-out, RMS in kcal/mol):

| Computed besides La | RMS |
|---|---|
| none (spread between ligands) | 1.68 |
| Dy | 0.28 |
| Sm or Eu, plus Tm | 0.15 |
| Nd, Dy and Lu | 0.12 (about experimental noise) |

Script: `series_rank.py`.

## Protein and peptide series fit the same two shapes (`protein_series_shapes.py`)

The data are from the curated table `workspaces/ln_selectivity_calibration_20261027/protein_ln_affinities.tsv`: 18 datasets, each with at least 6 lanthanides under one paper and one condition. Output is `workspaces/lady_series_20261028/protein_series_shapes.tsv`.

**Datasets that fit as well as NIST ligands do on the same metal subset.** "Shortcut" means reconstructing the rest of the series from La, Nd, Dy and Lu.

| Dataset | Two-shape residual (kcal/mol) | Shortcut residual (kcal/mol) |
|---|---|---|
| Mex-LanM, full 14-metal relative series (JACS Au 2024) | 0.12 | 0.17 |
| Mex-LanM (Cotruvo 2018) | 0.10 | — |
| Mex-LanD | 0.09 | — |
| LBT family (four datasets) | 0.14–0.21 | 0.14–0.30 |
| RTX | 0.08 | 0.12–0.21* |
| HEW5 | 0.11 | 0.12–0.21* |
| A0A7 | 0.18 | 0.12–0.21* |
| DFOB | 0.09 | 0.09 |

\* Range across RTX, HEW5 and A0A7.

**Misfits, both most likely measurement artefacts:**
- **GGR (Snyder 1990)** has a residual of 0.38. La–Gd fit to within 0.06 kcal/mol each. The heavy half zig-zags, with Tb +0.76 and Yb +0.79. Tb is the assay's luminescent reference ion.
- **Mex-LanM (Deblonde 2020, per-site value from β31)** has a residual of 0.43, with a heavy end that rises up to Tm. That contradicts both Cotruvo 2018 and the 2024 full series, which fall toward the heavy end and fit.

**Hans-LanM** has only La, Nd and Dy published (Mattocks 2023 SI Table S2), so the shapes cannot be tested on it.

## Hypotheses (fixed before any Nd/Lu job ran)

Each metal M is compared with La. EXP is NIST 46, from the same sources and ionic strengths as the La/Dy panel: acetate at I = 0.1, CDTA at I = 0.5.

- **S-H0, constant per-metal offset.** The offset c(M) comes from acetate. Pass: for each of Nd, Dy and Lu, the RMS residual over the six other ligands is below 0.8 kcal/mol.
- **S-H2, per-water bias.** The bias is c0(M) + b(M)·n_disp, calibrated on acetate plus EDTA. Pass: RMS over the other five ligands is below 0.8 for each metal. This mirrors H2 in `../lady_panel_20261028/README.md`.
- **Monotonic offsets.** This is a prediction of the aqua-leg explanation, not a usability gate: c(Nd) > c(Dy) > c(Lu), roughly in proportion to the radius contraction from La (0.053, 0.133 and 0.184 Å, Shannon CN9).
  - A non-monotonic pattern would mean the bias does not live in the aqua leg, even if offsets are ligand-independent.
- **Whole series.** Map the calibrated calculated values (Nd, Dy, Lu) to all 13 metals with the NIST basis, fitted without that ligand. Pass: RMS against NIST over all ligand × metal points below 0.8 kcal/mol, with acetate excluded because it is the calibration ligand.
  - The experimental targets differ in shape. Acetate, malonate and DTPA turn over (for example DTPA: Dy −4.56, Lu −4.05). EDTA, CDTA, TMDTA and NTA keep rising. So this tests the bend, not just the slope.

**Added 2026-10-29, before any Nd/Lu job ran** (H0 had already failed on La/Dy):
- **TTHA joins the ligands.** Like the panel, it is exempt from the CN test.
- **S-H3, stretched selectivity.** Per metal, calc = α(M)·exp + β(M), calibrated on acetate plus EDTA. Pass: RMS in experimental units over the other ligands below 0.8.
- **Scores exclude both calibration ligands** for S-H2, S-H3 and the whole-series test.
- **The whole-series test uses the form the panel's TTHA call picked** (from `PANEL_RESULT.json`): per-water → H2, stretch/constant → H3, since H0 already failed. The other forms are printed as well.
- **Secondary check (reported, not a gate).** The same forms, calibrated on the panel's four fit ligands (acetate, EDTA, NTA, DTPA) and scored on its four held-out ligands (malonate, CDTA, TMDTA, TTHA).
  - Why: for Nd, acetate and EDTA differ by only 1.2 kcal/mol experimentally (−0.41 vs −1.57), so the two-point H3 slope is poorly conditioned there.

## States

There are 33 Opt+Freq states. Each Nd state starts from its La analogue and each Lu state from its Dy analogue, with the metal swapped (`build_series.py`).
- **Geometry source:** the conformer-best structure for aqua, acetate and EDTA; the final panel geometry when it kept its CN; otherwise the panel's MACE start.
- **Nd:** aq8 and aq9; acmono w7 and w8; acbi w6, w7 and w8; edtax w2 and w3; and both hydration states for each panel ligand (DTPA w1 only).
- **Lu:** the same, except acetate mirrors Dy: acmono w7 and w8, acbi w6 and w7.
- **TTHA w0 for both** (added 2026-10-29). If its La/Dy analogue is unfinished, the start is the metal-specific MACE build in `ln_surrogate_20261027/nist/678/`.
- **Dropped 2026-10-29: Nd_tmdta_w3, Lu_cdta_w3, Lu_tmdta_w3 and Lu_mal_w7.**
  - La_tmdta_w3 lost its third water in DFT (CN 8). Dy ran only w2, w2 and w6 for these ligands.
  - So the smaller ion has no extra state worth testing. `build_series.py` moves any unrun dry-run copies into `workspaces/lady_series_20261028/dropped/`.

**Method:** the method block is byte-identical to the La/Dy runs. `basis_Nd.inc` and `basis_Lu.inc` use lcecp-1-TZVP with N_core 49 and 60, from `legacy/qmmm_lc/ecp_lib`, plus a pinned AutoAux J block. That block was generated with the same one-atom `NoIter PrintBasis` procedure; regenerating La with it reproduces the production `basis_La.inc` AuxJ exactly (130/130 lines).

## Launch (gated on the La/Dy panel)

Launch only if the La/Dy panel passes H0 or H2, or is borderline. If it fails both, the series run cannot be calibrated and is moot.

**Update 2026-10-29:** H0 has failed on La/Dy (see the panel README). The launch now depends on the frozen two-parameter forms (H2 per-water, H3 stretch): launch if at least one passes on the held-out ligands. The TTHA call picks the form at analysis time.

**Outcome 2026-10-29 19:15: not launched.** Neither frozen form passed. On CDTA and TMDTA the errors were H2 −2.98 and −2.46, and H3 −1.91 and −1.77 kcal/mol, against a bar of 0.8. The inputs and scripts below stay as a record. The proposed replacement is `../../LANTHANIDE_NIST_PREDICTOR_REVIEW_2026-10-29.md`.

The TTHA code changes are done (2026-10-29 13:00). They touch `build_series.py`, `series_analyze.py` and `submit_series.sh` (TTHA first in `ORDER`). Both scripts were dry-run in scratch: 33 inputs, and the analysis ran end to end on mock Nd/Lu states.

1. Rebuild inputs with the final panel geometries: `lnsel python build_series.py`.
2. Check queue load, then run `./submit_series.sh N_SE N_STD`, for example 4 2 while Jacob's queue is busy.
3. Run `workspaces/lady_series_20261028/watch_series.sh` in the background. It runs `series_analyze.py` at the end and writes `SERIES_RESULT.json`.
