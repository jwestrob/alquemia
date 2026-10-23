# Experimental target: intact, amidated LBT3

21 October 2026. Primary-source and repository evidence review only; no molecular calculation or job submission.

## Decision and implementation handoff

Use **NH3+–FIDTNNDGWIEGDELLA–CONH2** as the primary development construct. It is a 17-residue, one-site peptide; the terminal amide is not an eighteenth amino acid. With Asp3/Asp7/Asp13/Glu11/Glu14 deprotonated, its peptide charge is **−4** and either La(III) or Dy(III) complex has charge **−1**. This is the declared fixed protonation microstate, not proof that all protonation alternatives have negligible population.

The primary [9CEQ deposition](https://www.rcsb.org/structure/9CEQ) explicitly contains `FIDTNNDGWIEGDELLA(NH2)`, a terminal `NH2` residue, one La ion, and an author-assigned monomeric assembly. It supplies 20 solution-NMR conformers. The raw mmCIF and RCSB metadata are archived in the matching workspace. Existing 7CCO/7CCN structures remain useful starting conformations after explicit terminal-state preparation; 7CCN is Lu-conditioned, not an observed Dy structure. Different conformers are starting states, not independent validation cases.

Retain **NH3+–FIDTNNDGWIEGDELLA–COO−** as a separately named terminal-chemistry sensitivity: peptide −5, holo −2. It was the previous full-256 model's explicit assumption. Do not overwrite it, reinterpret its energies as amidated chemistry, or call it a proven preparation error.

## Separate experimental targets

### 2020 ITC: best verified thermodynamic development target

[Hatanaka et al., Scientific Reports 10, 19468](https://doi.org/10.1038/s41598-020-76527-y), Supplementary Table S1:

| Quantity | La | Dy |
|---|---:|---:|
| Binding ΔG, kcal/mol | −6.82 ± 0.03 | −9.18 ± 0.11 |
| Binding ΔH, kcal/mol | +3.83 ± 0.28 | +2.24 ± 0.12 |
| Fitted stoichiometry | 1.01 ± 0.05 | 0.97 ± 0.08 |

Conditions: **283.15 K, pH 6.0, 50 mM MES, 100 mM NaCl**; 15 µM peptide and 300 µM lanthanide-nitrate injectant. At least three experiments; the main figure labels errors as standard errors, but Table S1 does not independently specify the convention. Methods permit fixing stoichiometry at one for weak, nonsigmoidal curves. Apo/La/Lu gel filtration supports monomers; no separate Dy trace was verified.

Terminal chemistry is unresolved in this paper: Supplementary Fig. S5 assigns peptide charge −4, whereas the MD methods assign holo charge −2. No explicit synthesis terminus specification was recovered. Consequently the ITC result is a **conditional experimental benchmark**, not proof that either modeled terminal state is the assay construct.

### 2025 competition assay: explicit terminal chemistry, incomplete quantitative metadata

[KT et al., Chemistry – A European Journal 31, e202501318](https://doi.org/10.1002/chem.202501318), SI pp. S3, S6–S8, S23:

| WT quantity | La | Dy |
|---|---:|---:|
| Table S3 reported Kd/Ki, nM | 14202.6 | 111.3 |

Synthesis explicitly specifies free N-amine and C-terminal amide. The SI Materials sequence has **Y10**, conflicting with **I10** in the main text, computational methods and actual 9CEQ deposit. Treat Y10 as a likely textual error, with the residual assay-construct uncertainty recorded.

Tb luminescence competition used 50 nM peptide, 100 nM Tb, 50 mM MES/pH 6 and 100 mM NaCl; competitor IC50 values were converted using the reported Cheng–Prusoff relation. **Assay temperature is not explicitly specified**; 25 °C describes NMR and must not silently fill this gap. Methods mention 95% confidence intervals and Figure 6 describes standard deviations, but Table S3 provides neither numerically. The main-text La value is 14211 nM rather than Table S3's 14202.6; the table is the retained numeric source.

## What the simulation can be compared with

The balanced exchange convention remains:

`P–La + Dy(aq) ⇌ P–Dy + La(aq)`

`ΔG_exchange = ΔG_bound(La→Dy) − ΔG_aqueous(La→Dy)`.

Negative values favor Dy. For matched single-site dissociation constants this equals `RT ln[Kd(Dy)/Kd(La)]`. The ITC difference supplies a direct reported-free-energy comparison without inventing unreported Kd values.

Both measurements are conditional on solution composition and the peptide's accessible protonation/conformation ensemble. They are not affinities for one fixed donor geometry. Apo conformational contributions cancel in the balanced exchange when the same apo state ensemble and solution convention are used; different holo populations do not cancel. A binding simulation must integrate all relevant bound hydration and conformational states rather than retain the starting NMR restraints.

Explicit 0.1 M NaCl alone does not establish matching to a 50 mM MES assay. Buffer/anion complexation and metal activity corrections remain an observable-mapping task. If `alpha_M` relates unbound analytical metal to free aqua-ion activity, the simplified one-site relation is `Kd_conditional(M) = alpha_M * Kd_free(M)`, giving an exchange correction `RT ln(alpha_Dy/alpha_La)`. Whether that correction is negligible must be tested or bounded, not assumed. Protonation alternatives require their own state weighting.

The 2025 SI direct-binding equation prints a dimensionally inconsistent auxiliary expression (`c = Pt + x` where a product is expected in the quadratic). The actual fitting implementation/raw data were not recovered. Its competition model also needs checking at finite peptide/Tb concentrations. These are reasons to preserve assay-model uncertainty, not evidence that the published preference is wrong.

## Development, validation and exposure

- **Consumed development:** LBT3 2020 whole lanthanide table, 7CCO/7CCN, all historical LBT3 calculations, and 2025 WT affinity table. The 9CEQ structural evidence is new to this review but describes the same peptide, not an independent validation protein.
- **Consumed mechanistic comparison:** 2025 N5D. Its published qualitative loss of selectivity and Lu value were read in the main text. An attempted header-only SI extraction accidentally displayed the mutant Dy value and other non-La table rows. This is explicitly **not blind validation**. Mutant La was not displayed or transcribed, but withholding that single number does not restore blindness.
- **Future independent validation:** select and register another peptide/monomer before inspecting its affinity outcomes. A metadata-only review confirms that LBTv / 1TJB has a different amidated sequence and a monomeric assembly; see `FUTURE_VALIDATION_INVENTORY.md`. Its broad heavy-ion preference and published La/Tb values were already exposed by prior work, so it is a partially masked quantitative lead, not a fully blind holdout. A separate mutation panel whose La/Dy outcomes remain unopened is preferable if one can be verified from metadata. No independent validation case is falsely declared ready here.

## Progression status

**Ready:** explicit primary construct/state for preparation, preserved terminal sensitivity, exact consumed ITC values/conditions, raw primary sources and an exchange-observable contract.

**Still required before a quantitative experimental-validation claim:** terminal/assay mapping sensitivity, buffer/speciation treatment, usable uncertainty for any 2025 comparison, and a genuinely independent validation case. These do not prevent beginning interaction qualification or explicit-solvent implementation on the declared amidated state. They do prevent declaring a fully matched experimental prediction from one computed number.

No new structure prediction, simulation, fit, molecular calculation or scheduler submission was performed. Source downloads and text/metadata extraction are recorded in `EVIDENCE_INVENTORY.json`.
