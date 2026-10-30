# Constant-bias panel for the cluster-continuum La/Dy exchange (2026-10-28)

## Question

The careful cluster-continuum protocol (`../lady_relaxed_exchange_20261024`) over-favours Dy by the same amount on its two ligands:

| Ligand | Calculated | Experimental | Bias |
|---|---|---|---|
| acetate | −1.95 | −0.07 | −1.9 |
| EDTA | −6.35 | −4.01 | −2.3 |

Energies in kcal/mol. If the bias is ligand-independent (it would sit in the shared La(aq)/Dy(aq) leg), one experimental ligand calibrates it away and the protocol predicts relative La/Dy selectivity.

## Hypotheses (fixed before any panel result)

- **H0, constant bias.** Pass if both hold:
  - RMS of (calc − offset − exp) < 0.8 kcal/mol over the panel plus EDTA, with the offset taken from acetate;
  - the calculated order is CDTA < TMDTA < EDTA (experimental: −5.67, −4.87, −4.01).
- **H1, Born-like bias.** The bias tracks complex charge (offset + a·q_c²). That would be a two-parameter physics model.
- **H2, per-water bias** (added 2026-10-28 18:30, before any panel job finished). The bias grows with the number of first-shell waters the ligand displaces: offset + b·n_disp.
  - Motivation: `EXPERIMENT_THERMO.md` finds that the ΔH error grows with waters displaced (acetate −3.5, EDTA −8.6) and is partly cancelled in ΔG by TΔS.
  - In this panel, n_disp runs nearly parallel to ligand charge (acetate 1 → DTPA ~7.5). So H2 cannot be told apart from a bias linear in charge. It is distinguishable from H1, which goes as q_c².
  - If H2 wins the leave-one-out comparison, calibration needs two reference ligands, for example acetate and EDTA.
- **Neither.** Residuals erratic: static cluster-continuum is not usable at this precision.

## Ligands and states

The PBE0-D4/CPCM Opt+Freq method block is byte-identical to the acetate/EDTA runs. Aqua ions, water, acetate and EDTA are reused from the old runs.

| Ligand | NIST id | Complex charge | ΔΔG exp (kcal/mol) | States |
|---|---|---|---|---|
| malonate | 3480 | +1 | −0.655 | w6, w7 |
| NTA | 465 | 0 | −1.582 | w4, w5 |
| DTPA | 655 | −2 | −4.556 | w1 (q = 1 across the series) |
| CDTA | 600 | −1 | −5.674 (I = 0.5) | w2, w3 |
| TMDTA | 610 | −1 | −4.869 | w2, w3 |

**Starts** are the MACE-relaxed complexes from `ln_surrogate_20261027/nist` (`build_panel.py`). CDTA was checked to be trans.
- w3 = w2 plus one water, placed with the old `add_water`.
- Malonate w6 = w7 minus its most distant water.

**State choice:** the lowest-G water-balanced hydration state per metal, the same rule as `analyze.py`. A state whose coordination number changes during DFT is excluded.

**Limitation:** panel states are primary-only, with no conformer search. Acetate/EDTA are reported both conformer-best and primary-only. For acetate and EDTA, conformer search moved single states by up to 1.1 kcal/mol.

## Jobs and analysis

- **Jobs:** `workspaces/lady_panel_20261028/JOBS.txt` lists each state's Opt+Freq job and its chained afterany restart.
  - CDTA, DTPA and malonate run on standard_expand with 24 ranks.
  - TMDTA and NTA run on standard with 12 ranks.
  - `-N 1` throughout.
- **Watcher:** `watch_panel.sh` exits early if a running job's output is silent for 6 h, and otherwise runs `panel_analyze.py` once all jobs end.
- **Output:** `PANEL_RESULT.json`.

## Throttled schedule (2026-10-28, at the user's request: don't occupy the whole allocation)

The panel is capped at 4 standard_expand jobs plus 2 standard jobs running at once. Each slot is a dependency chain, and each step waits for the previous step's restart job.

| Slot | Chain |
|---|---|
| SE1 | La_cdta_w3 → La_mal_w7 |
| SE2 | Dy_cdta_w2 → Dy_mal_w6 |
| SE3 | La_dtpa_w1 → Dy_nta_w5 |
| SE4 | Dy_dtpa_w1 → La_nta_w5 → Dy_nta_w4 |
| STD-A | La_tmdta_w3, then nothing further on standard |
| STD-B | Dy_tmdta_w2, then nothing further on standard |

**Tier 1** is the literature-preferred hydration states above. For CDTA and TMDTA that means La q3 / Dy q2, as for EDTA.

**Tier 2** holds the alternate hydration states and runs only if tier 1 passes or is borderline:
- La_cdta_w2, Dy_cdta_w3
- La_tmdta_w2, Dy_tmdta_w3
- La_nta_w4
- La_mal_w6, Dy_mal_w7

Tier 2 was cancelled after about 20 minutes. Initial submission IDs are in `JOBS_initial_submission.txt`.

**2026-10-28 17:40: La_tmdta_w2 moved up from tier 2** (job 29036, restart 29037; standard, 12 ranks; a third standard slot).
- After 4 h and 24 TightOpt cycles, La_tmdta_w3's third water sat at 3.08 Å, with max gradient 4.5e-4. The other two waters were at 2.63 and 2.70 Å.
- That is next to the 3.1 Å CN cutoff. If the water leaves, the state is excluded. If it stays, a stretched water probably means w2 is the lower-G state.
- Without w2, the La-TMDTA result would rest on a state the rule may reject. The TMDTA leg of the order test would then be inconclusive.
- The other five tier-1 geometries keep their intended CN. At 4 h, the third water in La_cdta_w3 is at 2.78 Å.

## Interim result and TTHA addition (2026-10-29 11:30, before malonate, CDTA, TMDTA or TTHA finished)

**Finished ligands:** acetate, EDTA, NTA and DTPA. Calculated and experimental ΔΔG, in kcal/mol (conformer_best):

| Ligand | Calculated | Experimental | Bias | Waters displaced |
|---|---|---|---|---|
| acetate | −1.95 | −0.07 | −1.88 | 1 |
| NTA | −4.26 | −1.58 | −2.68 | 3.5 |
| EDTA | −6.35 | −4.01 | −2.34 | 6 |
| DTPA | −8.33 | −4.56 | −3.78 | 7.5 |

**H0 has failed, whatever the remaining ligands give.**
- After the acetate offset, EDTA, NTA and DTPA alone give RMS ≥ 0.86 over the six test ligands. The threshold is 0.8.
- The NTA residual (−0.80) can only become more negative. Dy_nta_w4 can only lower Dy's complex free energy.
- The calculated order of the four ligands matches experiment.

**Two post hoc forms fit the four ligands.** Both are frozen now and are scored only on held-out ligands (`FROZEN` in `panel_analyze.py`):
- **H2, per-water bias:** bias = −1.642 − 0.228·n_disp. In-sample RMS 0.41.
- **H3, stretched selectivity:** calc = 1.279·exp − 1.956. In-sample RMS 0.48.

**Malonate, CDTA and TMDTA can't separate the two forms.** Their predictions agree to within 0.5 kcal/mol. They still test whether a frozen two-parameter calibration transfers: pass if RMS < 0.8 in experimental units.

**TTHA (NIST 678) does separate them.** Jobs 29779/29780 (La) and 29781/29782 (Dy).
- TTHA is decadentate and leaves no water on either metal, so n_disp = 8.5, the most of any NIST ligand. Yet it is nearly unselective: experimental ΔΔG −0.68 kcal/mol.
- That combination is the regime of protein sites: many waters displaced, little selectivity.
- Predicted calculated values: H2 −4.26, H3 −2.83, H0 −2.56.
- **Call:**
  - TTHA calc ≤ −3.90 → per-water form.
  - TTHA calc ≥ −3.20 → stretch or constant form.
  - In between → inconclusive.
- These thresholds are for conformer_best. Primary-only intercepts are shifted by +0.141.
- Starts are the 9-donor, q = 0 MACE builds (`build_ttha.py`).
- TTHA has no water to lose, so the CN exclusion does not apply to it: 9 vs 10 bound donors are isomers of the same state.

**Why the form matters for the series and for proteins.** Every panel ligand prefers Dy, and almost every NIST ligand does: only 3 of 75 prefer La at Dy, and none by more than 0.1 kcal/mol. The two forms extrapolate differently to light-selective, many-water sites like the LanMs.

**Queue change:** La_mal_w7 (job 28736) was unchained from La_cdta_w3 so it could start at once. Jacob's other standard_expand jobs were user-held at the time.

## Dy_ttha_w0 restarted from La's geometry (2026-10-29 12:45, before any TTHA energy was final)

- **Problem:** the Dy MACE start had a strained TTHA backbone, with backbone angle gradients up to 0.14 Eh/rad (La's start: ≤ 0.012). TightOpt oscillated:
  - +0.26 Eh at cycle 2 and +0.04 Eh at cycle 4;
  - max gradient still 0.14 at cycle 7, where La's was 0.0017.
- **Same conformer as La:** heavy-atom RMSD 0.15–0.18 Å, with identical backbone dihedrals and donor set. So this was an optimizer problem, not a different state.
- **New start:** La_ttha_w0's TightOpt cycle-8 geometry with Dy swapped in. The method block is byte-identical. Jobs 30037/30038; the old run is in `workspaces/lady_panel_20261028/abandoned/Dy_ttha_w0_macestart/`.
- **Side benefit:** the TTHA pair now shares one conformer, so conformer noise largely drops out of the La/Dy difference. That matters because the call's inconclusive band is only 0.7 kcal/mol wide.

## Result (2026-10-29 19:15): the pre-registered test fails; the last three runs were cancelled

CDTA, TMDTA and Dy_nta_w4 finished (full output in `panel_result_20261029_1915.txt`). Dy_nta_w4 is now the lowest Dy NTA state, which moves NTA's calc to −4.60.

- **CDTA and TMDTA come out about twice as Dy-selective as experiment:** calc −11.66 and −10.45, against exp −5.67 and −4.87.
- **Frozen-form errors on them** (exp units, conformer_best), in the order CDTA, TMDTA:
  - H0: −4.11 and −3.70;
  - H2: −2.98 and −2.46;
  - H3: −1.91 and −1.77.
- **Every frozen form fails, whatever malonate and TTHA give.** With both of those at zero error, the best possible held-out RMS is 1.30 (H3). The order test (CDTA < TMDTA < EDTA) passes. The TTHA call between H2 and H3 is moot.
- **Cancelled at 19:20, restart job first, because no outcome could change the verdict:**
  - La_mal_w7 (28737/28736; Opt at cycle 214);
  - La_ttha_w0 (29780/29779; in Freq);
  - Dy_ttha_w0 (30038/30037; Opt cycle 35).

**Diagnosis from existing outputs** (`exchange_decomp.py`, `water_hbonds.py`). EDTA and CDTA have the same donors, charge, waters displaced and hydration pattern (La w3, Dy w2, n = 0). Experiment puts them 1.66 apart; the calculation puts them 5.31 apart. The exchange terms, in kcal/mol:

| | E | D4 | CPCM dielectric | rest of E | G − E |
|---|---|---|---|---|---|
| edtax | −6.20 | +0.51 | +16.74 | −23.45 | −0.02 |
| cdta | −10.72 | +0.57 | +5.35 | −16.63 | −0.80 |

- **The gap is not dispersion or thermal.** It is an 11.4 kcal/mol swing in the continuum term, only partly offset by the internal energy (+6.8).
- **Water H-bonds differ.** All three of La_cdta_w3's waters H-bond back onto carboxylate O (1.97–2.23 Å). La_edtax_w3 has one such H-bond and the Dy states have none.
- **Reading:** the exchange depends, by several kcal/mol, on how the explicit first-shell waters happen to be arranged in one optimized structure. A single-structure cluster-continuum calculation cannot average over that, and no ligand-level calibration can absorb it.
- **Post hoc, not a test:** a stretch refit on all six ligands (calc = 1.62·exp − 1.62) has a leave-one-out error of 0.79 kcal/mol in exp units. That is right on the bar, but only because the slope moved from the frozen 1.28. EDTA is the worst leave-one-out point (+1.29).
- **Untested alternative:** La_cdta_w3 is trapped in a high-energy water arrangement. A water-orientation search on La-CDTA (~1k CPU-h) would decide it; it has not been run.

## Next: whole-series extension (not launched)

Nd and Lu runs on this same panel are staged in `../lady_series_20261028/` (inputs, analysis and pre-registered criteria). They were gated on this panel passing, so they are not launched.

The proposed next step replaces the static protocol: a NIST-trained predictor with one computed ligand-strain descriptor. It is written up for outside review in `../../LANTHANIDE_NIST_PREDICTOR_REVIEW_2026-10-29.md` (2026-10-29; not started).
