# LBT3 exchange observable and assay mapping

21 October 2026. Scoped literature/mass-balance design; no molecular calls or submissions. The construct and exposure record in `REPORT.md` remain unchanged. No new validation affinities were inspected.

## Decision

**Calculate an explicitly defined La/Dy exchange free energy first, and report the conversion to the assay observable separately. MES cannot be assumed inert.** The present evidence can bound simple hydrolysis and finite-concentration assay effects, but cannot supply a defensible numerical La-versus-Dy MES correction under the exact assay conditions. Do not replace that missing correction with an invented confidence interval or fit it to the peptide result.

| Target | Verified solution/measurement | Meaning of retained result | Remaining mismatch |
|---|---|---|---|
| 2020 ITC | 283.15 K; pH 6; 50 mM MES; 100 mM NaCl; 15 µM peptide; 300 µM lanthanide nitrate injectant | Conditional binding ΔG difference Dy−La = −2.36 kcal/mol | Experimental terminal chemistry unresolved; total-metal titration does not establish free aqua-ion affinity |
| 2025 WT competition | pH 6; same nominal MES/NaCl; 50 nM peptide; 100 nM Tb; reported WT La 14202.6 nM and Dy 111.3 nM | Reported competition Ki values obtained with a fitted IC50 and Cheng–Prusoff conversion | Affinity-assay temperature and numerical confidence intervals absent; I10/Y10 synthesis-text discrepancy; free versus analytical Tb not specified |

The 2020 paper describes metal-into-buffer heat subtraction; that removes a calorimetric background, not automatically the equilibrium competition of MES with peptide. The 2025 direct-Tb fit uses an analytical-concentration quadratic, while competition uses an approximation. Neither source documents a conversion to free M3+ activity. These are therefore conditional assay targets. Primary methods: [Hatanaka 2020](https://doi.org/10.1038/s41598-020-76527-y) and [KT 2025](https://doi.org/10.1002/chem.202501318), archived with their SI in this workspace.

## 1. Define the species before comparing numbers

Exchange direction is `PLa + Dy(aq) -> PDy + La(aq)`; negative ΔG favors Dy. For the same peptide ensemble and solution standard state:

`ΔG_exchange,free = RT ln(Kd_free,Dy / Kd_free,La)`.

In a concentration convention at fixed medium, let `alpha_M = [unbound M, all soluble species]/[M3+]`. If only uncomplexed peptide and one metal-bound ensemble participate, with no ternary peptide–metal–buffer complex:

`Kd_cond,M = alpha_M Kd_free,M`

`ΔG_exchange,cond = ΔG_exchange,free + RT ln(alpha_Dy/alpha_La)`.

Thus stronger bulk complexation of Dy makes conditional peptide binding appear weaker; removing a positive bulk correction makes intrinsic Dy preference more negative. The aqua-ion state includes its equilibrium hydration shell, not a bare ion.

For a genuine thermodynamic activity convention, retain the activity-coefficient conversion separately. Equal +3 charge cancels leading charge-only terms, but does not establish cancellation of specific interactions. An explicit aqueous leg that already samples chloride association must not receive that chloride contribution twice. Record whether each simulation leg represents a constrained aqua species or the full unbound metal ensemble in its stated medium.

At fixed pH, peptide protonation and conformational states also form ensembles. The −4 amidated peptide/−1 complex is the declared initial microstate, not its measured population. Buffer effects on a peptide–metal–MES ternary state cannot be represented by bulk `alpha` alone.

## 2. MES: established concern, unresolved differential correction

[Mandal, Kretzschmar and Drobot 2022](https://doi.org/10.1007/s00775-022-01930-x) measured Eu–MES association using spectroscopy at room temperature and 150 mM NaCl. Their Table 1 gives **proton-release**, not ordinary deprotonated-ligand formation, constants:

| Reaction | Reported log10 constant |
|---|---:|
| Eu3+ + HMES -> EuMES2+ + H+ | −3.7 ± 0.1 |
| Eu3+ + 2 HMES -> Eu(MES)2+ + 2 H+ | −8.5 ± 0.2 |

Here the second product is Eu(MES)₂⁺ and the second constant is cumulative. The study used MES pKa 6.27. A concentration-based screening translation to pH 6/50 mM MES gives `HMES ≈ 0.0325 M` and `alpha_Eu ≈ 1 + 6.5 + 3.3 ≈ 11`. This illustrates substantial sequestration; it is **not** a measured La/Dy correction. The pH extrapolation, activity convention and temperature transfer need their own qualification. Preserve the source's uncertainty notation without relabeling it as a known confidence level.

[Azab et al. 2011](https://doi.org/10.1021/je200099n), primary Table 13, reports a 1:1 Dy–MES logK of 4.19 ± 0.02 at 25.0 ± 0.1 °C and 0.1 M KNO3; its table footnote labels uncertainties 3s. The same table gives Eu–MES 4.27 ± 0.03. Table content was inspected through an indexed copy of the primary PDF; a local original PDF was not retrieved. Protonation/formation conventions and different experimental media must be reconciled before using this value. A La–MES source exists, [Anwar and Azab 2001](https://doi.org/10.1021/je000187%2B), but only its primary abstract was accessible here; its numerical La constant is **not admitted**.

Consequences:

- No verified, matched La/Dy pair currently supports a best-estimate MES correction or a finite empirical bound on its difference.
- Include possible 1:2 complexes. Eu proves that restricting every metal to 1:1 is an additional assumption.
- Constants at 25 °C in nitrate medium are not exact 10 °C NaCl assay constants. Missing complexation enthalpies prevent a justified van't Hoff transfer.
- The often-quoted La value from a review was a search lead only; it is excluded from the executable input and all recommended corrections.

## 3. Hydrolysis and other solution chemistry

**Simple first hydrolysis is a small screening term at pH 6.** [Klungness and Byrne 2000](https://doi.org/10.1016/S0277-5387(99)00332-0), whose [institutional primary abstract](https://digitalcommons.usf.edu/msc_facpub/1723/) is archived, defines `beta1* = [MOH2+][H+]/[M3+]`. At 25 °C and zero ionic strength it reports logβ1* = −8.81 for La and −7.59 for Dy. Using these values at pH 6 gives hydrolysis-only alpha values approximately 1.0016 and 1.0257: a Dy−La exchange correction of only +0.014 kcal/mol at 298.15 K. This is a screening calculation, not the exact 10 °C assay. Their temperature measurements cover 25–55 °C; extrapolation to 10 °C must be labeled. Higher hydrolysis, precipitation and metal–peptide hydroxo states are not covered by this estimate.

**Chloride and nitrate require consistent bookkeeping.** Primary [Peppard et al. 1962](https://doi.org/10.1016/0022-1902(62)80109-2) directly measured weak chloride/nitrate associations, but at ionic strength 1.0 and without a Dy value in the retrieved series. [Bonal et al. 1998](https://doi.org/10.1039/A800428E) includes La and Dy nitrate calorimetry at ionic strength 2.0 mol/kg. Neither is an exact assay correction. NaCl is 100 mM, whereas nitrate delivered by the 2020 300 µM Ln(NO3)3 syringe is at most 0.9 mM before dilution. Acid used in stocks would change that bound; its amount was not recovered. Treat nitrate as a sensitivity until stock composition is verified, and avoid transferring high-ionic-strength constants without a medium correction.

**Other missing quantities:** air/CO2 exposure and dissolved inorganic carbon, trace phosphate/chelators, actual pH at measurement temperature, stock acidity and equilibration time. These may introduce complexes or precipitation; no contamination is asserted. MES contributes ionic species and counterions, so the ionic strength is not automatically exactly 0.100 M. State their absence or a controlled concentration in the simulation, rather than silently claiming assay identity.

## 4. The 2025 competition calculation can be checked without new chemistry

Use the conditional Kd measured in the same buffer for Tb. Do not apply a second MES correction to that already-conditional Tb constant.

For one competing site, let free peptide be `p`, analytical totals `Ptot`, `Ttot`, `Itot`, and conditional dissociation constants `Kt`, `Ki`. Exact soluble mass balances are:

`PT = p*Ttot/(Kt+p)`; `PI = p*Itot/(Ki+p)`; `Ptot = p + PT + PI`.

Solve the last equation for p and predict luminescence proportional to PT. This automatically includes peptide, Tb and competitor depletion. At zero competitor obtain `PT0`; at the exact half-signal point let `y=PT0/2`. Then:

`p = y*Kt/(Ttot-y)`; `A = Ptot-p-y`; `IC50 = A*(1 + Ki/p)`.

The published approximation instead uses `IC50 = Ki_reported*(1 + Ttot/Kt)` if its [Tb] means the added 100 nM. At 50 nM peptide/100 nM Tb/91 nM Kt, depletion is finite. The executable comparison below quantifies the resulting La/Dy ratio change, using both the methods' 91 nM and Table S3's 92.1 nM Tb constants. This is an idealized mapping sensitivity, **not a refit**: raw curves, backgrounds, fit implementation and numerical intervals are unavailable. The printed direct-binding auxiliary expression in the SI is dimensionally inconsistent; its actual software implementation remains unknown.

## 5. Required uncertainty and the next executable comparison

Keep three uncertainty classes separate: experimental repeat/fit uncertainty, simulation sampling uncertainty, and chemical/assay-model alternatives. The first two cannot absorb the third automatically. The reported 2020 La/Dy errors give about 0.11 kcal/mol in quadrature **if independent and with the same error convention**; they do not bound construct or buffer bias. No numerical 2025 confidence interval is invented.

At 283.15 K, an uncertainty factor of approximately 2.4 in `alpha_Dy/alpha_La` already means 0.5 kcal/mol; a factor of 5.9 means 1 kcal/mol. Changing −2.36 kcal/mol by its full magnitude requires a factor about 66 in the appropriate direction. These are decision thresholds, not claimed plausible ranges. One unknown differential proton uptake per exchange gives about 0.13 kcal/mol sensitivity to a 0.1 pH-unit shift at this temperature.

Prepared executable: `workspaces/lanm_selectivity_free_energy_20261021/experimental_target_v1/observable_mapping.py`, with `observable_mapping_inputs.json`. **Not executed here.** Root can run it as one CPU saved-data analysis on Slurm; no electronic/model calls or optimization are involved. It will report:

1. The retained 2020 conditional contrast and error-convention caveat.
2. 2025 exact mass-balance versus published-approximation sensitivity; 2025 temperatures are explicitly assumed scenarios, never recovered metadata.
3. Eu MES illustration and first-hydrolysis screening with source/convention labels.
4. A transparent alpha-ratio sensitivity grid and the correction needed to alter a proposed prediction by 0.5/1/2.36 kcal/mol.

**Next evidence needed for a matched prediction:** recover a compatible La/Dy MES speciation pair (including reaction definitions and temperature dependence), or measure a paired buffer-concentration series under the chosen assay. In the meantime, interaction qualification and explicit-solvent sampling can proceed for the declared intrinsic exchange model; report experimental agreement as conditional on the unresolved mapping. Do not reinterpret the present uncertainty as evidence that MES explains the entire mechanics failure.
