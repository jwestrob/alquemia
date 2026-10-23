# Observable mapping: actual results and independent collector review

21 October 2026. Root authorized local saved-data arithmetic with one thread. No molecular calls, scheduler submissions, new affinity outcomes or edits to the reviewed collector. Primary provenance and assumptions remain in `OBSERVABLE_CONTRACT.md` and `OBSERVABLE_EVIDENCE.json`.

## Result: finite-concentration competition is a small effect here

The idealized exact mass balances reproduce the reported IC50 values while changing the inferred conditional affinities as follows:

| Conditional Tb Kd used | Inferred La Ki, nM | Inferred Dy Ki, nM | Dy/La ratio |
|---|---:|---:|---:|
| 91 nM, SI methods | 13117.92 | 91.109 | 0.0069454 |
| 92.1 nM, SI Table S3 | 13111.36 | 90.991 | 0.0069399 |
| Published WT values | 14202.6 | 111.3 | 0.0078366 |

The resulting exchange shift is **−0.0679 to −0.0684 kcal/mol at assumed 283.15 K**, or **−0.0715 to −0.0720 kcal/mol at assumed 298.15 K**. Neither temperature is recovered metadata for the 2025 affinity experiment. This is a sensitivity calculation under one-site equilibrium, linear signal, and total-Tb interpretation of the published denominator; it is not a refit or a correction proven applicable to the original raw curves.

For the 91 nM Tb case, 22.928 nM of the initial 100 nM Tb is bound. At half signal, free peptide is 11.783 nM, peptide–Tb is 11.464 nM, peptide–competitor is 26.753 nM, and unbound Tb is 88.536 nM. An independent monotonic bisection solution reproduced all four cases; every concentration is positive and the largest checked mass-balance/equilibrium residual is **1.82×10⁻¹² nM**.

**Interpretation:** finite depletion is real but, under these assumptions, is too small to explain a multi-kcal/mol model discrepancy. The result does not resolve missing temperature, fitting uncertainty, terminal chemistry or the I10/Y10 assay-text discrepancy.

## Solution-speciation sensitivity remains separate

- The Eu–MES illustration gives alpha = **10.837**, with calculated aqua/one-MES/two-MES fractions **9.23% / 59.89% / 30.88%**. It uses the primary Eu proton-release constants and the stated concentration approximation. It is not an estimate of either La or Dy speciation.
- First-hydrolysis screening at 298.15 K and zero ionic strength gives alpha La **1.00155**, Dy **1.02570**, and a **+0.01412 kcal/mol** exchange correction. Exact assay conditions and higher species remain outside that screen.
- At 283.15 K, a factor **2.432 / 5.913 / 66.303** in alpha_Dy/alpha_La produces **+0.5 / +1.0 / +2.36 kcal/mol**, respectively. Reciprocal factors give negative shifts. The tested alpha-ratio grid is a sensitivity display, not an empirical uncertainty interval.
- The retained 2020 conditional contrast is **−2.36 kcal/mol**. Quadrature of its reported La/Dy errors is **0.114 kcal/mol only under the independence/same-error-convention assumption**; it excludes chemical-state and assay-model uncertainty.

**No La/Dy MES best estimate or finite empirical confidence bound is established.** Large absolute sequestration can cancel in a metal ratio, or it may not. Keep the intrinsic exchange and solution-to-assay conversion as separately reported quantities.

## Independent chemistry collector review: PASS for current artifacts

Reviewed `chemistry_comparison_v1/COLLECT.py` SHA `507f1f0c8968ea90fd164ef329763a59b8a32d8828ec3c9c5947a92b4706af10` and `pilot_partial_v2/RESULT.json`; all linked input pins were checked. Reconstructed the arithmetic independently without importing the collector.

1. **Cycle/sign:** donor contrast is `E(DyL)+E(LaW)−E(LaL)−E(DyW)`, the balanced reaction `LaL + DyW -> DyL + LaW`. Element and total-charge balance pass all 21 reported contrasts. Negative means Dy is favored in that donor environment relative to the specified water complex. For `aquo8_to_9`, negative means a more favorable eighth-to-ninth-water increment for Dy than La; this is a coordination increment, not a peptide preference.
2. **Units:** MACE eV and eV/Å are converted with 23.060547830619 kcal/mol per eV. ORCA Hartree energies use 627.509474; its Hartree/bohr gradients become forces via `−gradient × 627.509474 / 0.529177210903`. Tinker native kcal/mol/Å gradients are negated once. Saved reference conversions and Tinker signs replay exactly.
3. **Denominator:** force RMSE is over all **3N Cartesian components**, including metal and hydrogen atoms: acetate has 8 atoms/24 components; water has 4 atoms/12 components. The paired quantity is `(F_Dy−F_La)_model − (F_Dy−F_La)_reference` for each case. It is not metal-only error, RMS atom-vector magnitude, or a combined four-state exchange derivative. Maximum absolute error is also a component maximum.
4. **Reconstruction:** 21 contrasts agree within **1.02×10⁻¹⁰ kcal/mol**; 32 deformations replay exactly; all four force-error summaries agree within **3.56×10⁻¹⁵ kcal/mol/Å**. Paired model/reference La/Dy atom ordering and coordinates match.
5. **Coverage:** MACE 28, AMOEBA transfer 28, frozen-f PBE0 4, explicit-f PBE0 2. The explicit-f results are La-only, so no explicit-f La/Dy contrast or paired force error is available. Missing rows remain missing. The frozen-f central bidentate-acetate contrast is **−17.2616 kcal/mol**, versus MACE **−34.9017** and AMOEBA transfer **−31.1406**; these are fixed-configuration electronic contrasts with unresolved reference state/basis/accuracy limits.

No sign, unit or denominator defect was found in the current reported arithmetic. This review does not qualify the Hamiltonians, electronic states, sampling, or experimental selectivity.

## Artifacts and execution

Matching workspace `experimental_target_v1/` contains:

- `observable_mapping_result_v1.json`: actual arithmetic output, pinned source/input.
- `observable_mapping_validation_v1.json`: independent mass-balance checks and artifact pins.
- `chemistry_collector_independent_review_v1.json`: arithmetic, balance, unit and denominator checks with all reviewed input pins.

All work used one thread and zero molecular calls. Tool-observed mapping-process wall time was 0.124 s; collector-review process wall time was 0.130 s. CPU time was not separately measured. Earlier preparation records correctly remain marked unexecuted at their historical checkpoint; this report records the later authorized execution.
