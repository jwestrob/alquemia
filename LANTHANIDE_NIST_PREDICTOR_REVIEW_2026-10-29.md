# La–Lu complexation selectivity: what failed, and a proposed predictor (request for review, 2026-10-29)

**Request.** We want an outside critique before committing to the next step, and blunt answers are welcome. The most useful replies would point out:
- flaws in the proposed predictor or in its validation;
- an affordable physics route we have missed.

The specific questions are in section 7. This document is self-contained; section 8 gives the repository paths of the supporting files.

**Status:** this is a proposal only. Nothing has been run for it yet.

## 1. Goal

Predict, from ligand structure, how the stability of a 1:1 lanthanide complex changes across La→Lu in water at 25 °C. The method has to scale to thousands of ligands.

- **Quantity:** ΔΔG(Ln−La) = −RT ln10 · [log K1(Ln) − log K1(La)]. One log unit is 1.364 kcal/mol, and a negative value means the heavier ion is favoured.
- **Benchmark:** NIST SRD 46 v8, the Martell–Smith critically selected stability constants.
- **Worked example:** La/Dy.

## 2. The benchmark data

- **The series is close to one-dimensional.** For the 47 NIST ligands with complete La–Lu series, an uncentred SVD of ΔΔG(M−La) over 13 metals gives:
  - shape 1, a steady rise with the lanthanide contraction: 98.1% of variance;
  - shape 2, a bend at Gd/Tb: 1.3%;
  - shape 3: 0.3%.

  Knowing Dy−La reconstructs the other metals to 0.28 kcal/mol (leave-one-ligand-out), and knowing Nd, Dy and Lu gets to 0.12. Each ligand is therefore effectively one number, or two if the heavy-end bend matters.
- **Spread:** across 101 ligands, the sd of Dy−La is 1.12 kcal/mol.
- **Noise:** 33 ligands are listed at two or more ionic strengths.
  - Their per-value sd of log K(Dy) − log K(La) is 0.29 kcal/mol overall, and 0.57 for strong chelators (log K_La > 10).
  - This mixes real ionic-strength dependence with between-study scatter.

## 3. What has been tried

Every method was scored against NIST. None beats a simple donor-count regression.

| Method | Result |
|---|---|
| Gas-phase ML-potential (MACE-OMol) exchange energies for 75 ligands, used as a regression feature | Acts as a proxy for ligand charge. Adding it to donor counts makes leave-one-out (LOO) worse: 0.91 → 1.05 kcal/mol. |
| PBE0 + CPCM single points on the ML-relaxed structures | LOO 1.11 kcal/mol on its own, against 0.92 for donor counts + charge on the same ligands. Per-ligand scatter is ±3–4 kcal/mol. |
| Classical explicit-water alchemical FEP (12-6-4 ion parameters, OPC water, GAFF2-BCC; 13 λ windows × 6 ns) | Acetate −2.29 vs exp −0.07; EDTA −10.49 vs exp −4.01. La-EDTA unwraps: 71% of frames are tetradentate, though the crystal structure is hexadentate. |
| Static cluster-continuum DFT, tested with a pre-registered calibration panel (below) | Misses the held-out ligands by 1.8–4.1 kcal/mol. |
| Empirical ridge regression on 12 donor-group counts from SMILES, 101 ligands | LOO 0.85 kcal/mol, against 1.14 for predicting the mean. |

The force field and the static DFT are very different models, yet their error has the same sign per carboxylate. Both make carboxylate ligands too Dy-selective, and the excess (calc − exp) grows with the number of carboxylates:
- FEP: −2.2 kcal/mol for acetate, −6.5 for EDTA;
- static DFT: −1.9 for acetate, −2.3 for EDTA.

### The static DFT panel

**Protocol.**
- **Software and method:** ORCA 6.1.1; PBE0-D4/def2-TZVP; RIJCOSX; DefGrid3; tight SCF and optimisation; frequencies; CPCM water.
- **Metals:** large-core (4f-in-core) ECPs on La and Dy, so both are closed shell.
- **Exchange reaction:** La·L(H2O)b + Dy(H2O)c → Dy·L(H2O)a + La(H2O)d + n H2O.
- **State selection:** for each metal, the complex and the aqua ion take their lowest-G hydration state, with waters balanced. A state is dropped if its coordination number changed during optimisation.

**Calibration forms.**
- A constant offset fitted on acetate (H0) was registered before any panel result.
- H0 failed on the first four ligands (acetate, EDTA, NTA, DTPA). Two two-parameter forms were then fitted on those four and frozen before the held-out ligands finished:
  - offset + a per-displaced-water term;
  - a linear stretch, calc = 1.279·exp − 1.956.
- **Pass criterion:** held-out RMS below 0.8 kcal/mol, in experimental units.

**Result.**
- CDTA came out at −11.66 kcal/mol against −5.67 experimentally, and TMDTA at −10.45 against −4.87.
- Held-out errors were 1.8–4.1 kcal/mol for every form.
- Two ligands were unfinished (malonate, TTHA). Even if both were perfect, the best possible RMS is 1.30.
- The ordering CDTA < TMDTA < EDTA is right.

**Diagnosis.** EDTA and CDTA have the same donors, charge and hydration states (La with 3 waters, Dy with 2). Experiment puts them 1.66 kcal/mol apart; the calculation puts them 5.31 apart. Decomposing the exchange:
- Dispersion and thermal terms are the same for both.
- The difference is an 11.4 kcal/mol swing in the CPCM dielectric term, partly offset (+6.8) by the rest of the electronic energy.
- All three waters on La-CDTA H-bond back onto carboxylate O. La-EDTA has one such H-bond, and the Dy complexes have none.
- In a single optimised structure, the arrangement of a few explicit first-shell waters against the continuum moves the answer by several kcal/mol. In solution that arrangement is averaged, and no ligand-level calibration can absorb it.

**Two further points.**
- **Aqua ions:** separately, the cluster-continuum ΔΔG_hyd(Dy−La) is −59.4 kcal/mol, against about −67 experimentally.
- **Post hoc, not a pass:** refitting the stretch on all six ligands gives LOO 0.79, but only because the slope moved from 1.28 to 1.62.

## 4. Where the empirical baseline fails

**The error is concentrated.** The donor-count model already fits 94 of its 101 ligands to 0.54 kcal/mol. Seven ligands carry 63% of its squared error:

| Ligand | LOO residual (kcal/mol; + = more Dy-selective than predicted) |
|---|---|
| TTHA (10 donors) | −4.20 |
| CDTA | +3.21 |
| TMDTA | +2.35 |
| DTPA-bis(methylamide) | +2.23 |
| oxolane-2,3,4,5-tetracarboxylate | +2.00 |
| DTPA | +1.51 |
| EDTA | +1.43 |

All seven are polycarboxylates whose selectivity is set either by the backbone or by having more donors than a small ion can use. Donor counts encode neither.

**The clearest case is one donor set on different backbones.** NIST lists 15 diamine-N,N,N′,N′-tetraacetates with both La and Dy. Each has 2 N and 4 carboxylates; two carry one extra donor.

| Backbone | I (M) | log K La | log K Dy | dlogK |
|---|---|---|---|---|
| cis-1,2-cyclohexylene | 1.0 | 3.71 | 4.70 | 0.99 |
| trans-1,3-cyclohexylene | 1.0 | 7.02 | 8.41 | 1.39 |
| hexamethylene | 1.0 | 7.93 | 9.57 | 1.64 |
| trans-1,4-cyclohexylene | 1.0 | 7.30 | 9.03 | 1.73 |
| pentamethylene | 0.1 | 9.20 | 10.99 | 1.79 |
| cis-1,4-cyclohexylene | 1.0 | 7.13 | 8.96 | 1.83 |
| trans-1,2-cyclopentylene | 0.1 | 16.30 | 18.20 | 1.90 |
| cis-1,3-cyclohexylene | 1.0 | 6.24 | 8.18 | 1.94 |
| carboxyethylene (+1 carboxylate) | 0.1 | 14.30 | 16.27 | 1.97 |
| trans-1,3-cyclopentylene | 1.0 | 10.83 | 13.56 | 2.73 |
| ethylene (EDTA) | 0.1 | 15.36 | 18.30 | 2.94 |
| 2-hydroxytrimethylene (+1 OH) | 0.1 | 11.61 | 14.83 | 3.22 |
| 1,2-phenylene | 1.0 | 11.39 | 14.89 | 3.50 |
| trimethylene (TMDTA) | 0.1 | 11.22 | 14.79 | 3.57 |
| trans-1,2-cyclohexylene (CDTA) | 0.5 | 15.24 | 19.40 | 4.16 (3.68 at I = 1.0) |

- **The span is large.** The Dy−La values cover 4.3 kcal/mol, yet the donor model predicts 1.8–2.0 log units for all four of these that it contains.
- **Ring size matters.** Six-membered N–M–N chelate rings are more Dy-selective than five-membered ones (TMDTA and 2-hydroxytrimethylene against EDTA), as Hancock's ring-size rule predicts.
- **Long linkers cluster low.** Backbones too long to close a small ring between the two N cluster at 1.4–1.9.
- **The 2D graph misses geometry.** trans-1,2-cyclopentylene binds La about 1 log unit more strongly than EDTA does, but binds Dy no more strongly. The 2D graph cannot separate cases like these; the 3D backbone geometry can.
- **Most of the family is missing.** 11 of the 15 are absent from the current model: 9 because the PubChem name lookup failed, and 2 because they were measured at I = 1.0.

## 5. Proposal

Train a predictor on NIST for the one number per ligand. Its features encode chelate-ring size and backbone shape, plus one computed descriptor of 3D size fit. The data carry what we cannot compute reliably, mainly hydration, which varies smoothly with donor count. The computation carries what about 200 data points cannot teach: how a backbone fits a large ion versus a small one.

**Data: 101 → 208 ligands.**
- **58 La/Dy pairs whose structure lookup failed:** write the SMILES by hand, with stereochemistry, and check each against the NIST formula.
- **12 pairs measured at I > 0.5:** include them. Between ionic strengths, dlogK shifts by a median of only 0.09 log units.
- **37 ligands without La or Dy:** they have a light and a heavy Ln under the same conditions, so they enter through the series shape.

**Target.** The shape-1 amplitude, which is essentially Dy−La. Reconstructing the full series from it adds about 0.28 kcal/mol. If heavy-end selectivity matters, the shape-2 bend can be added as a second target.

**Features.**
- **From the 2D structure:**
  - donor-group counts and ligand charge;
  - donor excess beyond what one ion can bind;
  - counts of 5-, 6- and 7+-membered chelate rings closed by the bound donors;
  - backbone rigidity: rotatable bonds between donors, and whether the backbone is in a ring;
  - macrocycle cavity size.
- **One computed 3D descriptor, size-fit strain:**
  - Wrap the ligand around a dummy metal with metal–donor distances fixed at typical La values, relax it with a force field (MMFF), then repeat at typical Dy values.
  - Take strain(Dy) − strain(La), each measured from the free-ligand minimum. This takes seconds per ligand.
  - Lanthanide binding is ionic and non-directional, so distance constraints are a fair model of the metal.
  - No water and no continuum are involved, which are the components that failed above.
- **Excluded:** measured log K_La, because a new ligand will not have one.

**Model.** Ridge regression on about 20 named features, so each coefficient reads as kcal/mol per ring, per excess donor and so on. Nonlinear models serve only as a check.

**Validation.**
- **Baseline:** the donor-count model, refit on the same ligands and the same folds.
- **Headline number:** cross-validation that holds out whole ligand families, because plain leave-one-out lets near-identical analogues predict each other.
- **Block test:** all 15 diamine tetraacetates held out at once.

**Pass/fail, fixed before fitting.**
- **Pass:** family-held-out RMSE of 0.6 kcal/mol or less (about 70% of the variance), and 0.8 or less on the held-out block.
  - If the 2D features pass on their own, counting topology is enough.
  - If the strain term is needed, computed geometry carries the selectivity, and that is a physics component that scales.
- **Fail:** the data cannot teach backbone effects beyond counts and charge. The best available model stays near 0.85 kcal/mol. That separates strong selectors from weak ones but cannot rank close ones.

**Limits.**
- **Coverage is mostly carboxylates:** 84 of the current 101 ligands. Amides, phosphonates and macrocycles have fewer than a dozen examples each even after curation, and there are no hydroxamates. Predictions for those classes are extrapolation and get flagged by distance to the nearest training ligand.
- **Scope:** 1:1 complexes in water at 25 °C only.

**Cost.** About a day of work, mostly the structure curation. Compute is negligible: single-threaded, minutes.

## 6. Alternatives considered and set aside

- **Search water orientations for each static DFT structure.** This could fix CDTA for about 1k CPU-h. But it would have to be repeated for every complex, and it still freezes a single arrangement.
- **Explicit-water free energies with an ML potential.** For example, MACE-OMol (trained on ωB97M-V/def2-TZVPD) with La→Dy label-swap thermodynamic integration in a water box.
  - This is the one remaining route that combines sampled water with a DFT-quality metal–ligand interaction.
  - **Untried and unbenchmarked:** we have not tried it or measured its MD throughput.
  - **Sampling:** first-shell water exchange on these complexes takes nanoseconds or longer. The classical FEP could not sample EDTA's binding mode in 6 ns per window. We therefore expect weeks of GPU time per ligand.
  - **Risk:** it may inherit the per-carboxylate Dy overbinding that both other models showed.
  - **Cheap first gate:** a bare-aqua-ion test against ΔΔG_hyd(Dy−La) ≈ −67 kcal/mol.

## 7. Questions for reviewers

1. **Size-fit strain:** is the descriptor physically sound for Ln(III)? Is there a better-established way to compute size match, such as Hancock-style molecular-mechanics strain versus metal radius, or bite size? Which metal–donor distances and force-field settings would you use?
2. **Missing features:** which ones? Candidates include ligand basicity (Σ pKa, which would need predicted values for new ligands) and the change in coordination number across the series.
3. **Ionic strength:** is pooling ionic strengths acceptable at this noise level? Or should values be corrected (SIT, Davies) or restricted?
4. **Validation design:** with about 200 ligands, are family-held-out cross-validation and the block test the right design? Is 0.6 kcal/mol family-held-out RMSE the right bar?
5. **Prior work:** is there published QSPR or LFER work on lanthanide selectivity, as opposed to absolute log K, that we should benchmark against?
6. **Physics routes:** given the failure modes above, is there an affordable one we have missed, and at what cost? For example:
   - ML-potential MD with enhanced sampling;
   - QM/MM free energies;
   - a different continuum or cluster treatment;
   - a correction for the shared per-carboxylate error.
7. **Diagnosis:** do you disagree with any part of the failure diagnosis in section 3?

## 8. Where the evidence is (repository paths)

- `diagnostics/lady_panel_20261028/`: the panel README (pre-registration, result and decomposition), plus `panel_analyze.py`, `panel_result_20261029_1915.txt`, `exchange_decomp.py` and `water_hbonds.py`.
- `diagnostics/lady_explicit_solvent_20261026/RESULT_298K.md`: the classical FEP.
- `diagnostics/lady_relaxed_exchange_20261024/`: the static exchange protocol, the NIST SRD 46 decoder (`tools/nist46.py`) and the experimental thermodynamics (`EXPERIMENT_THERMO.md`).
- `diagnostics/ln_surrogate_20261027/`: the ML-potential and CPCM single-point regressions.
- `diagnostics/ln_selectivity_calibration_20261027/`: the NIST pipeline, the donor-count baseline and `DATA_AUDIT_20261029.txt`. The audit reproduces the noise estimate in section 2, every number in section 4, and the data counts in section 5.
- `diagnostics/lady_series_20261028/series_rank.py`: the series-shape numbers in section 2.
