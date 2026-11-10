# NIST charge response, La against Ca — 9 November 2026

**Result.** Three findings, in the order they were obtained.

1. **Frozen test: not supported.** The donor-count regression M1 predicts log K(La) − log K(Ca) for
   held-out ligand families with RMSE 1.58 log units against 1.97 for the training-fold mean
   (R² +0.13 [−0.62, +0.36]; n = 39 ligands, 9 families). M1 predicts log K(La) (R² 0.87) and
   log K(Ca) (R² 0.70) separately.
2. **A charge line exists across all pairs.** Over the 83 La/Ca pairs measured at one ionic strength,
   the difference rises 1.04 [0.69, 1.40] log units per unit of ligand negative charge
   (ρ +0.60 [0.41, 0.75]). The nitrogen-containing chelators carry it (ρ +0.84, n = 41). The 35
   oxygen-only ligands sit almost entirely at charge 1 or 2 and show ρ +0.10 [−0.27, +0.45].
3. **The line fails on proteins.** For eight protein rows in six families, the line's RMSE is 4.33 log
   units against 4.01 for the protein mean. Five families fall 1.9–7.3 log units below the line; the
   lanmodulins (coordination number 9–10) sit 1.2–3.6 above it.

Authorized by Jacob on 9 November (exploratory mining tracks). Finding 1 is the analysis subagent's
frozen design: `PLAN.md` sha256 `c2c04a1f…9246`, frozen 01:43 UTC 10 November before any fit.
Findings 2 and 3 are lead-agent additions outside that design (see DEVIATIONS.md). Written by the lead
agent from the subagent's reply and the saved tables; scores were checked against `summary.json` and
`addition_Mq.json`.

## Positive control

The historical pairing rule gives 171 La/Dy pairs (pairs.tsv: 171). The curation rule reproduces all
115 primary records. M1 on La/Dy scores 0.8545 kcal/mol (target 0.854).

## Pair counts

Ligands with a 1:1 constant for La and for the second metal at 25 °C (total / with feature vectors):
Ca 104/47, Sr 64/30, Mg 79/42, Ba 70/37, Mn 76/36, Zn 127/66, Cu 128/61, Cd 105/61, Al 43/25,
Fe³⁺ 56/28, Sc 21/13, Y 122/78, Lu 158/104. Feature vectors exist for the 115 La/Dy ligands only.
Scored set: same ionic strength, plus nearest-strength pairs within 0.1 M. Sr (23), Al (15) and
Fe³⁺ (9) fall below the frozen minimum of 25 and are unscored.

## Frozen scores

Family-held-out RMSE, log units. Mq = ligand net charge alone (added after the freeze).

| Target | n | M0 | Mq | M1 | M1 R² | RMSE(M0) − RMSE(M1), ligand interval |
|---|---|---|---|---|---|---|
| La − Ca | 39 | 1.97 | 1.95 | 1.58 | +0.13 | +0.40 [−0.02, +0.70] |
| log K(La) | 39 | 4.54 | 3.12 | 1.41 | +0.87 | +3.13, excludes 0 |
| log K(Ca) | 39 | 2.85 | 1.66 | 1.36 | +0.70 | +1.49, excludes 0 |
| La − Mg | 34 | 1.68 | 1.21 | 1.44 | −0.05 | includes 0 |
| La − Ba | 31 | 2.54 | 2.30 | 2.35 | −0.46 | +0.19 [0.03, 0.31] |
| La − Lu | 104 | 1.27 | 1.10 | 0.84 | +0.41 | +0.43, excludes 0 |
| La − Cd | 51 | 1.64 | 1.57 | 1.42 | −0.09 | excludes 0 |
| La − Y | 78 | 0.79 | 0.81 | 0.85 | −0.47 | M1 worse than the mean |

- Family intervals include 0 for every difference target.
- One-feature models on La − Ca: anionic oxygen donors 2.05, total donors 1.89.
- Sensitivity: same ionic strength only (n = 32) R² +0.11; all pairs (n = 47) R² −0.01.
- M1 coefficients for La − Ca (log units per group): amine N +0.75 [0.26, 1.23]; carboxylate
  +0.16 [−0.06, +0.40]; +0.35 per extra negative charge [0.08, 0.62].
- The held-out Mq residual follows the number of neutral donors (ρ +0.69 [0.48, 0.84]) and is flat in
  carboxylate count (ρ −0.01).
- La − Ca against the lattice-strain linear coefficient: ρ +0.54 [0.06, 0.84], n = 22.
- Penalty grid: interior for La − Ca; at the edge for the other targets; one decade wider leaves
  log K(La) and La − Lu unchanged.
- Noise floor: three ligands have pairs at two or more ionic strengths (pooled SD 0.35).

## Charge line on all pairs (lead addition)

Charge is read from the NIST protonation label (HnL neutral, so z = −n), which covers every ligand.
Ordinary least squares on q = −z; same ionic strength; 2,000 ligand bootstraps.

| Target | Ligands | n | Slope per charge | Intercept | ρ | Leave-one-out RMSE, line / mean |
|---|---|---|---|---|---|---|
| La − Ca | all | 83 | 1.04 [0.69, 1.40] | 0.61 | +0.60 [0.41, 0.75] | 1.64 / 2.19 |
| La − Ca | N-containing | 41 | 1.14 [0.76, 1.49] | 0.75 | +0.84 [0.72, 0.92] | 1.55 / 2.49 |
| La − Ca | O-only | 35 | 0.86 [−0.20, 2.01] | 0.84 | +0.10 [−0.27, 0.45] | 1.80 / 1.70 |
| La − Sr | all | 53 | 1.24 [0.87, 1.64] | 0.69 | +0.67 [0.45, 0.83] | 1.85 / 2.52 |
| La − Ba | all | 55 | 1.40 [0.92, 1.86] | 0.27 | +0.65 [0.39, 0.83] | 1.88 / 2.57 |
| La − Mg | all | 67 | 1.56 [1.04, 2.01] | −0.13 | +0.63 [0.44, 0.77] | 2.00 / 2.77 |
| La − Zn | O-only | 51 | 1.62 [0.79, 2.13] | −1.38 | +0.83 [0.74, 0.88] | 0.80 / 1.25 |
| La − Lu | all | 157 | −0.63 [−0.77, −0.48] | 0.16 | −0.59 [−0.70, −0.46] | 0.90 / 1.19 |
| La − Lu | O-only | 82 | +0.05 [−0.09, 0.21] | −0.61 | +0.08 [−0.12, 0.29] | 0.46 / 0.46 |

- **O-only coverage.** Charge 0: 1 ligand; 1: 15; 2: 16; 3: 1; 4: 2. Medians at charge 1 and 2 are
  2.28 and 1.96. The class mixes carboxylates with phenolates, hydroxamates, phosphates and hydroxide.
- **Class-held-out.** A line fitted on N-containing ligands predicts the O-only ligands with RMSE 1.54
  against 2.13 for the training mean; the reverse gives 1.62 against 2.78.
- **Other lanthanides against Ca** (all ligands, n = 36–45): slopes run from 1.08 (Nd) to 1.77 (Yb);
  Pr 1.34, Eu 1.29, Tb 1.44, Dy 1.63, Er 1.68, Lu 1.35.
- **Within the series.** Among O-only ligands, La − Lu is −0.61 at every charge (slope +0.05). The
  charge dependence of the La-to-Lu slope lives in the N-containing chelators.
- **Zn and Cd.** Over all ligands La − Zn and La − Cd are flat in charge (ρ −0.01, −0.10), because
  the highly charged ligands carry amine nitrogens.

## Protein comparison (lead addition, post hoc)

Protein pairs come from `~/Documents/Berkeley/Nikasha/protein_ln_ca_table_2026-11-09/` (Mac). q is the
number of acidic residues in the first shell. Lines: La 0.61 + 1.04 q; Pr 0.31 + 1.34 q;
Eu 1.05 + 1.29 q.

| Protein | Ln | q | Anionic O | CN | Observed | Line | Residual |
|---|---|---|---|---|---|---|---|
| Glucose/galactose receptor (2FW0, Ca) | La | 3 | 4 | 7.0 | −1.46 | 3.73 | −5.20 |
| α-Lactalbumin (1F6S, Ca) | La | 3 | 3 | 6.8 | −3.54 | 3.73 | −7.27 |
| Calmodulin (1PRW, Ca) | Eu | 3.5 | 4.5 | 7.0 | +0.26 | 5.56 | −5.31 |
| RTX domain (5CVW, Ca) | La | 2 | 2.9 | 6.75 | +0.80 | 2.69 | −1.88 |
| LanP (9VY7, La) | La | 4.8 | 6.2 | 8.8 | +1.10 | 5.60 | −4.50 |
| Hans-LanM (8DQ2, La) | La | 4 | 8 | 10.0 | +5.95 | 4.77 | +1.18 |
| Mex-LanM pH 5 (8FNS, Nd) | Pr | 3.75 | 6 | 9.25 | +7.30 | 5.34 | +1.96 |
| Mex-LanM pH 7.2 (8FNS, Nd) | La | 3.75 | 6 | 9.25 | +8.13 | 4.51 | +3.62 |

- RMSE of the line 4.33; RMSE of the protein mean 4.01; mean residual −2.17; Spearman between line
  and observation +0.35 (p 0.40).
- Coordination number against the observed difference: ρ +0.77 (nominal p 0.025, n = 8 rows, two of
  them one protein). Coordination number against log K(Ln): ρ +0.86.

## Extended protein comparison (added 10 November UTC, post hoc)

A literature gap-fill (`~/Documents/Berkeley/Nikasha/protein_ln_ca_gapfill_2026-11-09/`, Mac; every
new value a single-agent extraction, most from abstracts) adds 16 matched pairs, four of them with a
structure: troponin C high- and low-affinity sites (Eu, +2.37 and +2.42), carp parvalbumin (La, +1.08)
and bovine trypsin (Gd, +1.20).

- Twelve structured rows in seven families: RMSE of the line 3.90 against 3.30 for the protein mean;
  mean residual −2.38; Spearman +0.27 (p 0.39). All nine rows outside lanmodulin fall below the line,
  by 1.9 to 7.3 log units.
- The eight rows with a Ca-bound structure (coordination number 6–7) span −3.54 to +2.42, median
  +0.94. The first table's four such rows (−3.5 to +0.8) were the low end of that range.
- All 28 pairs: lanmodulin +5.95 to +8.13 (three rows); the 25 pairs in 11 other families have a
  median of +1.10 (range −3.54 to +3.20).

## Interpretation (inference)

In flexible chelators, charge rises together with denticity and chelate-ring count, and every added
arm can wrap either ion. In the extended set, the eight rows whose sites hold Ca at coordination
number 6–7 show differences of −3.5 to +2.4 (median +0.9) at two to four acidic residues, where the
line expects +2.7 to +4.9. The lanmodulins supply 9–10 donors and exceed the line. LanP limits a pure
capacity reading: coordination number 8.8 with 6.2 anionic oxygens and a difference of +1.1. What
separates lanmodulin from the other 11 families is therefore unexplained by charge, donor count or
coordination number in these data.

For mining, the supported statement is about affinity: deliverable donor count (anionic oxygens,
total coordination) ranks log K(Ln) at ρ +0.86 in this set. La-over-Ca selectivity remains a separate
term that neither ligand charge nor a small-molecule scale predicts for proteins.

## Limits

- The O-only ligands barely sample charge 3 and above, so NIST cannot calibrate sites with four to
  six carboxylates.
- Donor classes for the all-pairs line come from name keywords: 92% agreement with the feature table
  on 100 ligands; seven of eight disagreements are nitroso-naphthols.
- Six of seven nearest-strength pairs put La at 0.1 M against Ca at 0 M, consistent with a −0.4 to
  −0.9 log bias (inferred); the same-strength set gives the same reading.
- Protein rows mix evidence classes, three are cross-study, multi-site proteins use site averages,
  and the comparison was designed after the protein table was read.
- Unverified: original-study provenance of each constant; the species behind the large DTPA and EDDS
  values.

## Reserved material left unopened

`workspaces/ln_structure_increment_20261030/exclusions.tsv`; the SpyCI-LAMBS and PLM paths; every
`lit_untrusted/` directory; the protein and peptide series tables.

## Files

- `diagnostics/nist_charge_response_20261109/`: `PLAN.md`, `extract_pairs.py`,
  `fit_charge_response.py`, `add_charge_baseline.py`, `charge_line_all_pairs.py` (sha256
  `48c7a3b4…a87bc`).
- `workspaces/nist_charge_response_20261109/` (untracked): `la_pairs_all_metals.tsv`,
  `pair_counts.tsv`, `la_ca_all_conditions.tsv`, `y1_coefficients.tsv`, `y1_primary_predictions.tsv`,
  `y1_Mq_residuals.tsv`, `summary.json`, `addition_Mq.json`, `charge_lines.tsv`,
  `charge_line_ligands.tsv`, logs.
