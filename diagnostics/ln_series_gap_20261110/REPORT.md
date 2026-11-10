# The lanthanide series splits into two branches at the half-filled shell — 10 November 2026

**Result.** In NIST log K1 series, each ligand needs a cubic polynomial in ionic radius
(log K = c0 + c1·r + c2·r² + c3·r³, four numbers per ligand), and all ligands then share one
correction to the size axis: an extra gap between Gd and Tb of 31 mA on the Shannon nine-coordinate scale
(95% interval 18.5 to 61.5; about 2.2 element steps). The series is two smooth branches, La–Gd and Tb–Lu. The
gap was found after the frozen readings were scored, so it was tested on 24 ligands that no fit had touched:
with the gap fixed at 31 mA in advance, leave-one-metal-out error fell from 0.209 to 0.134 kcal/mol
(difference 0.075, interval 0.050 to 0.099; better for 21 of 24 ligands), and those ligands alone give a gap of
29.5 mA (16.5 to 76.5) at the same position.

Three further findings:

1. **Protein and peptide series show no gap.** Twelve datasets pooled: +2 mA (−18 to +8.5). The one 14-metal
   lanmodulin series: +5 mA (−3.5 to +16.5). Small-molecule corrections should stay out of protein predictions.
2. **Yttrium takes the heavy branch only with oxygen donors.** With oxygen-only ligands it sits with Dy and Ho,
   its size twins (fraction of the gap 0.18, −0.10 to 0.52; replicated on the confirmation ligands, −0.22,
   −0.52 to 0.04). With nitrogen-containing ligands it sits on the La–Gd side (1.14, 0.68 to 2.50; 14 ligands).
3. **The break is a step, located after Gd.** A step after Gd beats the spin-pairing cusp at Gd of f-electron
   theory as a one-number shape, on both ligand sets, and beats a step one position either side.

Yesterday's "ligand-independent per-metal pattern of about 0.1 kcal/mol" is this structure seen through a
quadratic polynomial: most of it is the average of ligand-specific r³ terms, and the rest is the Gd|Tb gap.

"Quadratic", "cubic" and "quartic" below always mean polynomials of degree 2, 3 and 4 in the size axis.

Self-chosen by the lead agent (Claude Fable) under Jacob's instruction of 10 November ("work on whatever you
like"). Mac, CPU only, about three minutes of compute. `PLAN.md` (sha256 `c8ae855e…3b08`) was frozen before any
fit; `ADDENDUM_confirmation.md` (sha256 `364c5f06…4f50`) before the 24 confirmation series were read. The
literature subagent's reply supplied the theory tables and the Kawabe values. Figure:
`figures/lanthanide_two_branches.png`.

## Frozen readings of PLAN.md (47 complete-series ligands)

Score: leave one ligand out for everything shared, then leave one metal out for the ligand's own numbers
(kcal/mol; bootstrap over ligands).

| Model | Shared numbers | Own numbers per ligand | Score |
|---|---|---|---|
| Q, quadratic in Shannon radius | 0 | 3 | 0.204 [0.170, 0.239] |
| O, quadratic + shared offset | 11 | 3 | 0.176 [0.149, 0.201] |
| C, quadratic in a shared coordinate | 12 | 3 | 0.190 [0.163, 0.217] |
| OC, both | 23 | 3 | 0.178 [0.153, 0.204] |
| S, two free shapes | 22 | 3 | 0.173 [0.146, 0.200] |
| K, cubic in Shannon radius | 0 | 4 | 0.158 [0.124, 0.191] |
| Q4, quartic | 0 | 5 | 0.176 [0.136, 0.217] |
| T, quadratic + two f-electron shapes | 0 | 5 | 0.185 [0.146, 0.225] |

| Reading | Outcome |
|---|---|
| 1. Shared fine structure exists | **Pass.** Q − O = +0.028 [+0.009, +0.045] |
| 2. Offset or coordinate | **Undetermined.** O − C = −0.014 [−0.036, +0.005] |
| 3. One coordinate suffices | **Fails** for C against S (+0.017 [+0.009, +0.025]); **passes** for OC against SO (+0.007 [+0.000, +0.013]) |
| 4. Both terms | **No.** OC does not beat O (−0.003 [−0.017, +0.009]) |
| 5. f-electron shapes against generic terms | **Undetermined.** Q4 − T = −0.009 [−0.038, +0.021]. Share of the mean residual in the theory plane 0.66 (exact p = 0.008); a cubic alone explains 0.86 |
| 6. Amplitudes rise with N donors | **Pass as frozen** (A1 ρ = +0.60 [+0.36, +0.77]; A3 +0.40 [+0.14, +0.61]); the trend vanishes given denticity (partial +0.13 [−0.17, +0.44]) |
| 7a. Yttrium has one place on the coordinate | **Fails.** A single position (0.33–0.35) does not beat the Dy twin (0.32) |
| 7b. Yttrium's position rises with N donors | **Pass.** ρ = +0.62 [+0.32, +0.82] (Q), +0.68 [+0.40, +0.85] (C); same with N atoms |
| 8. Transfer to 18 protein datasets | **Fails.** O better in 11 of 18 (sign test p = 0.48); C in 8 of 18 |

The per-ligand cubic, which the plan carried only as a comparator, beat every shared model. That redirected
the analysis.

## Post hoc: the gap (47 ligands)

- Cubic in a free shared coordinate (12 numbers, ends pinned): 0.136 [0.104, 0.167], against 0.158 for the plain
  cubic (difference +0.023 [+0.006, +0.040]) and 0.138 for three free shapes (33 shared numbers).
- A shared offset adds nothing to a cubic (K − KO = +0.001 [−0.001, +0.004]); the coordinate beats the offset
  (+0.022 [+0.004, +0.039]).
- Fitted positions minus Shannon (mA): La 0, Ce +3.8, Pr +7.5, Nd +10.7, Sm +12.8, Eu +13.0, Gd +16.2,
  Tb −8.3, Dy −11.4, Ho −10.5, Er −8.2, Tm −9.6, Yb −5.9, Lu 0. Displacement relative to the two neighbours:
  Gd +14.3 [+5.6, +20.7], Tb −10.7 [−21.6, −4.6]; every other metal within ±2.6.
- One number in place of twelve: light branch +31.0 mA relative to the heavy branch [+18.5, +61.5]. Score 0.138
  [0.109, 0.168]; gain over the plain cubic +0.020 [+0.009, +0.032] (interior metals only: +0.007
  [+0.002, +0.014]). Gd|Tb is the best single break position in 945 of 1,000 resamples.
- Residual removed by one break at each neighbour pair (%): Gd|Tb 12.8, Tb|Dy 6.3, Yb|Lu 5.4, all others ≤ 2.4.
- Other base scales: Shannon CN8 +28.9 mA (2.1 element steps [1.3, 3.4]); atomic number 1.2 steps [0.5, 2.2];
  Shannon CN6 (values recalled, unverified) 0.9 steps [0.2, 12], on a scale whose own Gd–Tb step is already
  15 mA against 12–13.
- Shape of the feature, each as one shared number on a cubic per ligand:

| Shape | Residual SS (none 6.141; free coordinate 4.873) | Score, 47 | Score, 24 (number fixed from the 47) |
|---|---|---|---|
| Step after Gd | 5.358 | 0.138 | 0.134 |
| Step after Tb | 5.753 | 0.145 | 0.157 |
| Step after Eu | 6.140 | 0.158 | 0.207 |
| Spin-pairing tent, cusp at Gd (E1) | 5.621 | 0.149 | 0.155 |
| Tent with Kawabe's (q + 25) scaling | 5.640 | 0.149 | 0.156 |
| Triangle cusp at Gd | 5.954 | 0.153 | 0.129 |
| E3 shape | 6.051 | 0.157 | 0.220 |

  Tent minus step: +0.010 [+0.002, +0.019] on the 47 and +0.020 [+0.012, +0.029] on the 24.

## Confirmation on 24 untouched ligands (frozen addendum)

Twenty-two of the 24 lack Ce; 20 are oxygen-only; 16 list Y.

| Test | Outcome |
|---|---|
| C1. Gap fixed at +31 mA | **Confirmed.** K 0.209, KS 0.134; K − KS = +0.075 [+0.050, +0.099]; better for 21 of 24. Interior metals only: +0.014 [+0.008, +0.020] |
| C2. Free coordinate from the 47 | **Confirmed.** K − KC = +0.083 [+0.053, +0.110]; KC beats KS by 0.008 [+0.002, +0.012] |
| C3. Independent estimate | **Replicates.** Gap +29.5 mA [+16.5, +76.5]; Gd|Tb is the best break position (1,594 of 2,000 resamples) |
| C4. Offset adds nothing to a cubic | **Not replicated.** K − KO = +0.016 [+0.012, +0.019]; a fifth of the gap's gain |
| C5. Yttrium | **Replicates for oxygen-only ligands** (heavy branch better: −0.065 [−0.084, −0.042]; best fraction −0.22 [−0.52, +0.04]; n = 15). One nitrogen-containing ligand, consistent with the light side; untestable |
| C6. Proteins, 12 datasets with ≥ 9 lanthanides | **No transfer.** KS better in 5 of 12 (p = 0.77); with a quadratic the gap makes 10 of 12 worse (p = 0.04) |

On these 24 ligands the quadratic with the 47-ligand offset scores 0.128, level with the free coordinate
(0.127). Their slopes are small and similar, so an additive pattern and a slope-proportional one cannot be told
apart there. The ligand-class difference in gap size suggested by the 47 (oxygen-only 18 mA, nitrogen-containing
43 mA) is not supported: the oxygen-only confirmation ligands give 62 mA [41.5, 82].

## Yttrium (42 + 16 ligands)

- Model-free: with oxygen-only ligands yttrium's closest lanthanides are Dy and Ho (RMS difference 0.17
  kcal/mol each; Tb 0.18). With nitrogen-containing ligands the closest is Tb (0.40; Gd 0.49, Dy 0.50, Ho 0.63).
- Against a lanthanide of its own radius on each ligand's smooth curve, yttrium is weaker by +0.08 kcal/mol
  (SE 0.03) with oxygen-only ligands and by +0.46 (SE 0.09) with nitrogen-containing ligands (cubic baseline);
  by N atoms: 0 → +0.08, 1 → +0.22, 2 → +0.60. Spearman with N atoms +0.66 [+0.38, +0.82].
- Predicting log K(Y) from each ligand's lanthanides alone: on the light side of the gap 0.296; on the heavy
  branch 0.368; difference +0.072 [+0.012, +0.140]. By class the answer reverses, as in finding 2 above.

## Proteins and peptides

Own gap of each dataset with an F-based 95% interval (cubic): lanmodulin on agarose, 14 metals, +5.0 mA
[−3.5, +16.5]; every other dataset has an interval wider than 80 mA. Pooled over the 12 datasets: +2.0 mA
[−18.0, +8.5] (bootstrap over datasets); without the six LBT datasets +2.5 [−11.5, +8.5]. For comparison, the
71 NIST ligands one at a time: median +36 mA; the interval excludes zero on the positive side for 20 and on the
negative side for 6.

Yttrium in proteins: two public series list it. In the *E. coli* galactose-binding protein it lies +0.16
kcal/mol from a lanthanide of its radius and closest to Ho. In the 14-metal lanmodulin series it is +0.67
kcal/mol weaker than its size twin and closest to Tm and Er: lanmodulin disfavours yttrium by more than its
size explains.

The protein table also holds a nine-metal lanthanide series for a thrombospondin-repeat calcium-binding domain
(Chem. Sci. 2025, doi 10.1039/d5sc02315g).

## Candidates tested for the origin, and what each must explain

- **Kawabe's hydration stabilisation of the light aqua ions** (1999, the published ligand-independent
  correction): cosine with the observed shared pattern +0.32 (p = 0.31); adding it lowers no score and worsens
  the cubic (K − KH = −0.021 [−0.035, −0.007]). Not supported by these data.
- **Spin-pairing (E1) cusp at Gd:** captures two thirds of the gap's gain and loses to the step on both sets.
- **A change of coordination number after Gd** gives a continuous free energy with a change of slope; the fitted
  coordinate shows opposite displacements at Gd and Tb, the signature of a step.

Any explanation has to give a step between f7 and f8 that scales with each ligand's size selectivity, that
f-electron-free yttrium shares with oxygen donors, and that protein sites lack. I have none that does all four
(inference: the protein result points at the small complexes themselves, for example their inner-sphere water,
and away from the aqua ion they share with proteins).

## What this changes in practice

1. **Small-molecule numbers overstate light-versus-heavy selectivity for proteins.** The project's small-molecule
   target, ΔΔG(Dy−La), straddles the gap. Share of a chelator's fitted selectivity that comes from the gap
   (ligands with |ΔΔG| ≥ 0.5 kcal/mol; `scripts/gap_share.py`): La→Dy median 25% (IQR 21–31%, 44 ligands);
   Nd→Dy 30% (28–32%, 20); La→Lu 24% (19–32%, 49). Protein series show no gap, so a surrogate or series shape
   fitted on NIST carries a term proteins lack. Fill in a protein's series with a smooth polynomial in Shannon
   radius.
2. **Scoring size-only physics methods.** A method in which the lanthanides differ only by size gives a smooth
   series and misses the step by construction: about 0.1 kcal/mol for acetate-class ligands and about 1 kcal/mol
   for EDTA-class. Score small-molecule benchmarks within a branch, or add the gap term.
3. **Protein REE profiles** should compress to three or four polynomial coefficients plus one yttrium offset.
   A clearly negative fitted gap means the assay's readout is relative to a competing chelator.
4. **Inference, untested:** in protein-based separations the light/heavy cut at Gd|Tb comes from the chelator
   eluent; proteins have to supply selectivity within a branch.

## Deviations from PLAN.md

1. **N-donor feature.** `aminoacid_N + amine_N` counts carboxymethyl arms and double-counts nitrogen. Readings 6
   and 7b are reported as frozen and with N atoms from the molecular formula; the outcomes agree.
2. **Everything under "Post hoc"** (cubic families, the gap, shape comparison, yttrium branches, ligand classes,
   protein gaps, other base scales) was added after the frozen readings were scored. The gap, the free
   coordinate, the break position and the oxygen-only yttrium result were then tested on untouched ligands
   under `ADDENDUM_confirmation.md`.
3. **Yttrium in models O and OC.** The first run computed apparent positions from slopes of mean-centred series;
   those two lines are void. Prediction errors (reading 7a) are unaffected.
4. **External size scales.** The D'Angelo 2011 radii and the Persson 2008 distances are paywalled; the planned
   comparison of the coordinate with aqua-ion sizes was not made.
5. **Theory shapes.** The tables match Kawabe 1992 Table 1. Kawabe scales each term by (q + 25); the scaled
   shapes were added as a sensitivity run (T*: 0.170).
6. **Bootstrap of the coordinate** used two starts (Shannon and the full-data optimum); 1,000 resamples.
7. **Frozen prediction file** for the sealed ortholog profiles: `FROZEN_PREDICTION.md` (P1: no gap; P2: yttrium
   as a size twin). A dry run on public protein series, made after hashing, supports P1 (−1.0 mA, −15 to +9) and
   puts P2 in doubt for lanmodulin (+0.67 kcal/mol in the one public series); the note is appended to the file.

## Limits

- 71 ligands, all carboxylates or aminocarboxylates, one database, mostly ionic strength 0.1 M.
- The gap's interval is wide on the upper side (18.5 to 61.5 mA).
- The protein conclusion rests on one precise dataset and eleven imprecise ones.
- The nitrogen-containing yttrium result has 14 ligands and no independent replication.
- Shannon CN6 radii were typed from memory.
- The physical origin of the step is unexplained.

## Files

- `PLAN.md`, `ADDENDUM_confirmation.md`, `PLAN.sha256`. `FROZEN_PREDICTION.md` (hashed in `PLAN.sha256`) stays with the private tables.
- `scripts/`: `lib.py`, `run_main.py` (frozen readings), `run_proteins.py`, `run_followup.py`,
  `fit_cubic_coordinate.py`, `fit_step.py`, `explore_groups.py`, `yttrium_branch.py`, `run_confirm.py`,
  `robustness.py`, `gap_share.py`, `make_figure.py`, `sealed_gap_test.py`.
- `results/` (NIST-derived; kept outside this repository, with `data/` and `figures/`): `run.log`, `summary.json`, `followup.log`, `confirm.log`,
  `robustness.log`, `coordinate.tsv`, `coordinate_cubic.tsv`, `per_ligand.tsv`, `per_ligand_followup.tsv`,
  `protein_gap.tsv`, `step_fit.json`, and the matching summaries.
- `figures/lanthanide_two_branches.png` and `.pdf`.
- `data/` with `SHA256SUMS`; `prior_scripts/` (copies of the 28 October and 9 November scripts).
