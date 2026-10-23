# Where the fixed-La force discrepancy occurs

**The discrepancy includes systematic water attraction and donor-specific external force/torque errors. It persists after internal molecular force distortions are excluded.** This is disagreement with the consumed uncorrected ECP-PBE0 surface, not proof that the published La model is implemented incorrectly or that PBE0 is physical truth. The Dy-only pilot is held; no new energy or force calculation was performed.

Saved-only analysis retained all14 La fixtures and every covalent body, replayed both stage RMS values and all category bookkeeping, using one CPU thread for0.108s. Result: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/matched_dry_v1/la_diagnosis_v1/RESULT.json`, SHA256 `9dd4768a672e13e66b36451fb131b8256fe485dd18da7011833e23cd8ab195bb`. It includes every original vector, body membership, radial axis, tangential residual and torque.

## Full fixture ranking

Category-balanced rigid RMS in kcal/mol/Å, with torque divided by1Å:

| Fixture | RMS | Dominant visible pattern |
| --- | ---: | --- |
| Matched dry mono-acetate | 14.557 | Ligand torque40.998kcal/mol and substantial nonradial force |
| Hydrated mono-acetate −0.15Å | 12.238 | Ion force, ligand torque33.951 and opposing water/ligand radial residuals |
| Hydrated acetamide −0.15Å | 10.741 | Ligand inward radial residual−36.601 |
| Hydrated mono-acetate central | 9.832 | Ligand torque23.704, ion radial residual−21.165 |
| Matched dry acetamide | 9.020 | Near-radial ligand residual−17.902 |
| Matched dry bi-acetate | 8.961 | Near-radial ligand residual−17.888 |
| Water dimer | 8.153 | Water inward radial residual−16.306 |
| Hydrated mono-acetate +0.15Å | 7.764 | Ligand torque15.250 and ion radial residual−18.327 |
| Hydrated acetamide central | 6.744 | Ligand inward radial residual−24.394 |
| Hydrated acetamide +0.15Å | 4.870 | Ligand inward radial residual−15.109 |
| Hydrated bi-acetate −0.15Å | 4.562 | Ligand radial residual−11.668 plus water inward response |
| Hydrated bi-acetate central | 4.393 | Ligand radial residual−12.146 plus water inward response |
| Hydrated bi-acetate +0.15Å | 4.388 | Ligand radial residual−13.176 plus water inward response |
| Aquo8 | 3.083 | Nearly balanced ion force hides water radial residuals−10.552 to−10.668 |

Residual means AMOEBA minus reference. For ligand/water, positive radial points away from the ion; negative means more inward attraction or less outward repulsion than the reference. The table's RMS mixes the explicitly declared translation/rotation scales; it is not an atomic-force RMSE or a thermodynamic error.

## Systematic versus donor-specific response

**Water:** every represented water has a negative outward radial residual, including the dimer, aquo8 and all hydrated donor/radial fixtures. Aquo8 water net-force errors are about10.57–10.68kcal/mol/Å even though the ion residual norm is only about0.533, because opposing water directions largely cancel at the ion. In hydrated fixtures the water radial residuals span roughly−3.82 to−19.38. Thus a metal-force-only fit would miss a systematic water response error. Internal water bending cannot explain this net-body force by itself.

**Amide:** the ligand residual is predominantly radial and inward, from−36.60 through−24.39 to−15.11 over the hydrated radial sequence; matched dry is−17.90. At the hydrated center its Cartesian net-force residual is (−20.608,0.000013,−13.093)kcal/mol/Å. This is an external ligand response, with a smaller but nonzero torque residual7.842kcal/mol.

**Bi-acetate:** hydrated ligand inward residuals remain about−11.67 to−13.18 across the radial sequence. Matched dry is−17.89. Central Cartesian ligand residual is (−7.945,−7.646,−5.701)kcal/mol/Å. Tangential and torque residuals are smaller than in mono-acetate.

**Mono-acetate:** the direction changes with hydration. The hydrated ligand residual points outward (+10.59,+8.06,+4.89 across−/0/+), while the matched dry ligand is inward−12.67 with tangential norm11.26. Its torque residual is large:40.998 dry and33.951/23.704/15.250 hydrated. At the hydrated center the ligand vector is (+9.332,+0.141,−0.033), while the ion vector is (−21.364,+1.292,−5.327). This cannot be represented as a single uniform radial attraction correction; environment and orientation matter.

## Decision limits

The independent source audit found the expected published AMOEBA22 La block, SI water/donor typing and pair corrections; these saved arrays do not identify an implementation error. Published La calibration used CP-DLPNO-CCSD(T)/ZORA and hydration evidence, whereas this comparison uses uncorrected frozen-f ECP-PBE0-D4. Reference-method sensitivity and the known published mono-acetate limitation remain plausible causes. Diagnose those before releasing a Dy-only fit. The coarse acetamide-La derivative flag stays unchanged, and no conclusion about analytic derivative failure follows from that one0.15Å secant. Reserved protein/orientation validation outcomes remain unopened.
