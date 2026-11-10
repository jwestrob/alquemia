# One coordinate for the lanthanide series? — PLAN (frozen 2026-11-10, before any fit)

Self-chosen by the lead agent (Claude Fable) under Jacob's instruction of 10 November ("work on whatever you
like"). Mac only, CPU only, no cluster jobs. NIST-derived tables stay in this private folder.

## Question

The lattice-strain test of 9 November (`ln_lattice_strain_20261109`) left one loose end: after a per-ligand
quadratic in Shannon radius, every ligand shares a per-metal residual pattern of about 0.1 kcal/mol whose origin
was unverified. Two readings of a shared pattern make different, testable predictions:

- **Offset.** A free-energy term that belongs to the metal's aqua ion and is added to every ligand's series
  unchanged: y = quadratic_L(r_M) + c_M.
- **Coordinate.** The radius scale itself is slightly off: there is one coordinate ρ_M such that every series is a
  quadratic in it, y = quadratic_L(ρ_M). The residual each ligand shows then scales with that ligand's local slope
  and reverses sign on the far side of a ligand's optimum.

If the coordinate reading holds, the 47 ligand series collapse onto parabolas in one data-derived coordinate, and
that coordinate can be compared with independent measurements of aqua-ion size. Yttrium, which has no f electrons
and was left out of every earlier fit, then gets a position on the same coordinate from one or two shared numbers.

## What I already know (stated so the plan is honest about its blindness)

From the 9 November report: with the fit forced through La, the mean residual is Ce −0.10, Pr −0.12, Nd −0.10,
Sm −0.07, Eu −0.03, Gd +0.09, Tb +0.07, Dy +0.08, Ho +0.09, Er +0.08, Tm −0.01, Yb −0.09, Lu −0.10 kcal/mol; a
cubic explains 0.90 of it and a step at Gd 0.60; adding it as a shared offset lowers the three-metal prediction
error from 0.186 to 0.165. I have seen no fit with a free intercept, no coordinate fit, no yttrium value beyond
the two table rows printed while locating the files, and none of the external radius tables requested from the
literature subagent (I recall the D'Angelo 2011 aqua radii only roughly).

## Data

- `data/surrogate_nist.tsv` (sha256 in `data/SHA256SUMS`), parsed as `prior_scripts/series_rank.py`: 47 ligands
  with complete La…Lu series (14 metals; Pm absent), log K1 from NIST SRD 46. y_LM = −1.3643 · log K1 (kcal/mol,
  298 K), with a free intercept per ligand in every model. 42 of the 47 also list Y.
- Secondary set: the 71 ligands with at least 10 of the 14 metals (shared parameters fitted with missing cells).
- Features from the same table: N donors = `aminoacid_N` + `amine_N` + `pyridine_N`; `carboxyl`; `q_lig`;
  denticity `nlig_La_cn8`.
- Protein and peptide series: `data/protein_ln_affinities.tsv`, datasets built as
  `prior_scripts/protein_series_shapes.py` (18 datasets). Sealed tables are not touched.
- Shannon CN9 radii (Å): La 1.216, Ce 1.196, Pr 1.179, Nd 1.163, Sm 1.132, Eu 1.120, Gd 1.107, Tb 1.095,
  Dy 1.083, Ho 1.072, Er 1.062, Tm 1.052, Yb 1.042, Lu 1.032.

## Models

Every model gives each ligand three free numbers (intercept and two loadings) unless noted. H(x) is the projector
onto span{1, x, x²}; ȳ is the mean series over ligands.

| Model | Form | Shared degrees of freedom |
|---|---|---|
| Q | quadratic in Shannon CN9 radius | 0 |
| O (offset) | Q + c_M, c = (I − H(r)) ȳ | 11 |
| C (coordinate) | quadratic in ρ; ρ_La and ρ_Lu pinned to Shannon; ρ minimises Σ_L ‖(I − H(ρ)) y_L‖² | 12 |
| OC | quadratic in ρ plus c_M (C fitted to series centred on ȳ) | 23 |
| S (rank-2 ceiling) | intercept + two free shapes (top two principal directions of the metal-centred series) | 22 |
| SO | same about ȳ | 33 |
| K | cubic in Shannon radius, 4 numbers per ligand | 0 |
| Q4 | quartic in Shannon radius, 5 per ligand | 0 |
| T (f-electron shapes) | Q + A1_L·t1 + A3_L·t3, 5 per ligand | 0 |

C is nested between Q and S: it asks that the two shapes be a vector and its square.

**Fitting C and OC.** Minimise tr[(I − H(ρ)) Σ_L y_L y_Lᵀ] over the 12 free ρ (scipy, analytic profile over ligand
parameters). Starts: Shannon CN9; Shannon CN8 rescaled to the pinned ends; linear in atomic number; first
principal direction rescaled; 20 jittered copies of the CN9 start (σ = 0.01 Å, seed 42). Report the number of
distinct optima and whether ρ stays monotonic. No monotonicity constraint.

**Theory shapes for T** (Jørgensen's refined spin-pairing energy theory; coefficients of the Racah parameters E¹
and E³ for the 4f^q ground terms, q = 0…14; Pm dropped):
t1 = (9/13) · [0, 0, 1, 3, 6, 10, 15, 21, 15, 10, 6, 3, 1, 0, 0]; t3 = [0, 0, 9, 21, 21, 9, 0, 0, 0, 9, 21, 21, 9, 0, 0].
Both are positive, so a positive amplitude means the complex is destabilised relative to the aqua ion where the
aqua ion gains most from spin pairing (Gd) or from orbital correlation (Nd, Ho, Er): the direction a more covalent
ligand should give. The literature subagent is checking these tables against a primary source; a correction would
be a recorded deviation.

## Scores (kcal/mol)

- **S_LOMO (primary).** Leave one ligand out to fit everything shared; for the held-out ligand leave one metal out,
  fit its own numbers on the other 13, predict the 14th. RMS over 47 × 14 predictions. Applies to every model.
- **S_4 (practical).** Shared parts from the other 46 ligands; the ligand's three numbers from La, Nd, Dy, Lu;
  RMS over the other 10 metals. Three-number models only.
- **S_in.** Shared parts from the other 46; the ligand's numbers from all 14 metals; RMS residual.
- Uncertainty: bootstrap over ligands (2,000 resamples, `numpy.random.default_rng(42)`), paired for differences,
  95% percentile intervals. Shared fits are computed once per held-out ligand.
- Robustness (reported, no reading): each ligand's series divided by max(SD across metals, 0.5 kcal/mol) before
  the shared fits; and the 71-ligand set.

## Frozen readings

Differences are in S_LOMO unless stated.

1. **Shared fine structure exists** if the best of O, C, OC beats Q with the interval of the difference excluding 0.
2. **Offset or coordinate.** S_LOMO(O) − S_LOMO(C): interval above 0 reads "coordinate"; below 0 reads "offset";
   otherwise "undetermined".
3. **One coordinate suffices** if the upper end of the interval of S_LOMO(C) − S_LOMO(S) is at most +0.015.
   The same statement for OC against SO.
4. **Both terms** if OC beats both O and C with intervals excluding 0; otherwise the simpler winner stands.
5. **f-electron shapes.** S_LOMO(Q4) − S_LOMO(T): interval above 0 reads "the theory shapes beat generic smooth
   terms of equal count"; below 0 reads the reverse; otherwise "undetermined". Also reported: the share of the
   mean Q residual (free intercept) lying in the plane of (I − H) t1 and (I − H) t3, with its exact p value
   (a random direction in the 11-dimensional residual space puts a Beta(1, 4.5) share there).
6. **Covalency direction.** Spearman correlation of A1_L (and of A3_L) with the number of N donors over the 47
   ligands; predicted positive; bootstrap interval; "supported" if the interval excludes 0 on the positive side.
7. **Yttrium.** Fit one shared position ρ_Y (models C, Q; in O and OC also one shared c_Y) on the other ligands
   that list Y, with each ligand's quadratic fixed from its 14 lanthanides, and predict Y for the held-out ligand.
   (a) "Yttrium has a place on the coordinate" if that prediction beats both fixed twins (Y = Ho and Y = Dy) with
   intervals excluding 0. (b) Apparent position per ligand by one Newton step from the shared position,
   ρ_Y(L) = ρ_Y + (y_LY − ŷ_LY)/slope_L, for ligands with |slope_L| ≥ 5 kcal/mol/Å at ρ_Y. Classical claim
   (Siekierski; Bau 1996): yttrium behaves like a lighter, larger lanthanide the more covalent the ligand.
   Prediction: Spearman correlation of ρ_Y(L) with N donors is positive; "supported" if the interval excludes 0.
8. **Transfer to proteins.** For each of the 18 protein and peptide datasets, with every shared quantity taken
   from the NIST fit and nothing shared refitted: leave-one-metal-out RMS under Q and under the winning shared
   model. "Transfers" if the median paired difference is positive and a two-sided sign test gives p < 0.05.

## Descriptive, no gate

- The fitted coordinate with bootstrap intervals, as a deviation from Shannon CN9; its comparison with Shannon
  CN8 and CN6 radii and with aqua-ion size scales from the literature subagent (D'Angelo 2011 radii; EXAFS mean
  Ln–O distances), each rescaled to the pinned ends. Statistic: cosine between deviation vectors over the 12 free
  metals, with the exact random-direction p value (cos² ~ Beta(1/2, 11/2)).
- Per metal, across ligands: regression of the Q residual on the ligand's local slope (intercept = offset-like
  part, slope = coordinate-like part).
- Vertex and curvature of each ligand in the fitted coordinate; a data-collapse figure.
- A frozen prediction file for the sealed ortholog profiles, written last and hashed, for Jacob to run or ignore.

## Outputs

`scripts/` (analysis code), `results/` (tables, `run.log`, `summary.json`), `figures/`, `REPORT.md`,
`DEVIATIONS.md` if needed. Code, plan and report can go to the public repository later; tables derived from NIST
values stay here.
