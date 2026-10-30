# Ln surrogate: enthalpy/entropy split and series structure (2026-10-27)

**Option 1: MACE gas dE_cn8 vs ddH and TddS (n = 20 overlap)**
- dE tracks ddG (r 0.78), but only through charge (charge vs ddG r −0.87).
- Partial r(dE, ddH | charge) = +0.20; partial r(dE, TddS | charge) = +0.13.
- **Verdict: dead.** Gas MACE dE is a charge proxy and carries no bonding signal beyond charge.
- ddH and TddS are uncorrelated across ligands (r 0.08).
- CDTA vs EDTA differs mainly in ddH (−0.4 vs +1.7): preorganisation acts through enthalpy.

**Option 2: lanthanide-series structure (NIST, 105 ligands with a complete La..Lu series)**
- PC1 explains 91.6% of the variance, PC2 6.9%. PC2 is the early-vs-late (Gd break) shape.
- Donor counts predict PC1 (LOO R² 0.57) but not PC2 (−0.14).
- Predicting Dy−La from a single measured X−La (LOO, NIST): Gd 0.39 kcal, Eu 0.42, Lu 0.48, Nd 0.64. Baseline sd 1.15.
- **Transfer to proteins**, using the line fitted on NIST and applied unchanged to the protein table:

  | Probe | RMSE (kcal) | Observed sd (kcal) | n |
  |---|---|---|---|
  | Tb | 0.28 | 1.26 | 13 |
  | Gd | 0.27 | 1.42 | 10 |
  | Eu | 0.38 | 1.21 | 8 |
  | Nd | 1.57 | 2.03 | 17 |

  Nd fails on the LanMs, whose preference reverses late in the series (conformational).
- **Caveat:** Tb and Gd sit next to Dy, so this is near-interpolation. It is not a structure-based prediction.

## Full MACE CN8 fit (76/76 built)
- Leave-one-out RMSE (kcal): donor counts 0.91, charge 0.96, dE + charge 0.98, donors + dE + charge 1.05. sd 1.19.
- **Dead.**

## MACE series-shape scan (La, Nd, Gd, Dy, Ho, Lu; jobs 27955–27960)
- **Kill criterion (stated in advance):** the partial r of the curvature given charge must exceed 0.3.

| Quantity | n | r | Partial r given charge | p |
|---|---|---|---|---|
| Curvature (Gd − Lu/2) | 74 | +0.33 | +0.26 | 0.023 |
| Slope (Lu − La) | 74 | +0.53 | +0.28 | 0.016 |
| Curvature, CN-conserved ligands only | 54 | | +0.26 | 0.054 |

- **Verdict: fails, narrowly.** There is weak signal, and the calculated spread is 8× the experimental one.
- Gas relaxations drop waters for the smaller ions: 19 of 75 ligands lose coordination number by Lu.
- LOO RMSE for Dy−La: donors + charge 0.93 → 1.04 with the calculated terms added.
- LOO RMSE for curvature: 0.42 → 0.37.

## Next: CPCM PBE0 large-core single points on the MACE La/Dy CN8 structures (jobs 28106–28113)

## Hydration literature (agents; notes in workspaces/ln_selectivity_calibration_20261027/lit_untrusted/)

**Ion hydration entropy**
- ΔhydS(Dy³⁺) − ΔhydS(La³⁺) = −68 J/K/mol (Rizkalla & Choppin 1991, ±40–50%).
- At 298 K that is TΔΔS ≈ +4.8 kcal/mol at full dehydration.

**Water count on the metal in the complex (q)**
- EDTA: La–Nd q = 3; Dy–Lu q = 2.
- CDTA: switches earlier than EDTA (Gd is 82% q = 2).
- DTPA, DTPA-BMA: q = 1 across the series.

**Dehydration step** (Janicki & Mondry 2019): [Gd(EDTA)(H2O)3]⁻ → q2 + water.
- ΔH = +18.2 kJ = +4.35 kcal; TΔS = +66 J/K × 298 K = +4.7 kcal. So ΔG ≈ −0.35 kcal.
- Mean q along La→Lu tracks the measured complexation entropy with r = −0.99 (EDTA and CDTA).

**Implication**
- The large EDTA/CDTA ΔΔH and TΔΔS mostly reflect the q3→q2 equilibrium.
- That equilibrium is near-thermoneutral in ΔG; it is centred mid-series by construction.
- EDTA and DTPA have similar ΔΔG (−4.0 vs −4.6 kcal) but TΔΔS +5.7 vs +2.4 kcal.
- So ΔΔG lives mainly in the fixed-q exchange term, not in the entropy.
- This corrects the earlier framing that selectivity is driven by charge-neutralisation entropy.

## CPCM incident

| Problem | Cause | Fix |
|---|---|---|
| 7 of 8 pools exited after one task | ORCA's mpirun consumed the pool's stdin task list | ORCA now reads stdin from /dev/null |
| Pool 28111 hung 8.5 h in the RI/J step | Job spanned 2 nodes | `-N 1` added |

Resubmitted as jobs 28536–28543.

## CPCM PBE0 large-core on MACE CN8 La/Dy structures (152/152, jobs 28536–28543)

**Leave-one-out RMSE (kcal; n = 75; sd 1.19)**

| Model | RMSE |
|---|---|
| donors + charge | 0.92 |
| dE_cpcm | 1.11 |
| dE_cpcm + charge | 0.94 |
| all | 1.00 |

**Verdict: FAIL** by the pre-stated criterion.

**What improved over gas MACE**
- Magnitudes are now physical: EDTA −2.3 vs −4.0 exp; DTPA −4.7 vs −4.6.
- Raw order CDTA < EDTA < TTHA is correct (−3.1, −2.3, +2.6), but TTHA has the wrong sign.

**What remains wrong**
- Per-ligand scatter is ±3–4 kcal: acetate +4.0 vs −0.1; NTA +0.6 vs −1.6.
- DOTA −22 (structure failure).
- Within charge class, r = 0.09 (q1, n = 32), 0.62 (q2, n = 27), 0.14 (q4, n = 10).
- There is a systematic offset of about +1.6–2.3 kcal at q1/q2.
- r with ddH is 0.42 (n = 20), so it fails the ΔΔH criterion as well.

## Reassessment (2026-10-28): what each result actually tested

**Careful cluster-continuum** (lady_relaxed_exchange_20261024: PBE0-D4/CPCM Opt+Freq, crystal/CREST starts, each metal in its preferred CN)

| Ligand | Calculated | Experimental | Bias |
|---|---|---|---|
| acetate | −1.95 | −0.07 | −1.9 |
| EDTA | −6.35 / −5.89 | −4.01 | −2.3 / −1.9 |

- The bias is ~constant for a 1-donor and a 6-donor ligand.
- Calibrating on acetate predicts EDTA at −4.0 to −4.5.

**Same nominal exchange** (mono-acetate, CN8 throughout), PBE0/CPCM single points:
- −1.0 on solvent-optimised geometries (ingredient tests);
- +4.0 on MACE gas geometries (NIST run).

So the NIST CPCM failure was a geometry failure.

**State pairing alone** moves the acetate exchange from −1.0 to −12.0 (ingredient tests).

**Cluster-continuum aqua leg:** ΔΔG_hyd(Dy−La) = −59.4 (La aq9, Dy aq8; per-water standard state RT ln 24.46 + RT ln 55.34), against −67 to −69 experimental.
- The model under-hydrates Dy³⁺, plausibly through the missing size-dependent second shell.
- Most of this cancels in exchanges; about 2 kcal does not.

**12-6-4 FEP:** a force-field result (metal–carboxylate terms transferred from water; La-EDTA 71% tetradentate). It is not a test of explicit-solvent physics.

**Retracted:** "selectivity is driven by charge-neutralisation entropy that the models exclude".
- EDTA's large TΔΔS is the q3→q2 switch, with ΔG ≈ 0.
- Explicit solvent contains solvent entropy.

**Next decisive test:** is the bias constant across ligands? Same careful protocol on:

| Ligand | Complex charge |
|---|---|
| malonate | +1 |
| NTA | 0 |
| DTPA | −2 |
| CDTA | −1 |
| TMDTA | −1 |

**Pass:** residual RMS < 0.8 kcal/mol after the acetate-calibrated offset, and CDTA < TMDTA < EDTA ordered correctly.

**Housekeeping:** jobs 25164/25165/25202 were hung for 51 h (multi-node ORCA) and were cancelled. All ORCA scripts now carry `-N 1`.
