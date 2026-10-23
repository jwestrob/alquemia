# Fixed metal-polarizability response: no usable correction

## Decision

**Close this finite polarizability probe without another fit.** All eight actual probes retain failed scientific checks. The largest change in the W8→W9 paired energy error is only 0.0332 kcal/mol, compared with the 4.3383 kcal/mol baseline error. The extra local response directions do not produce a physically admissible correction in the unconstrained projection.

This does not rule out every distant nonlinear parameter combination. It does rule out treating the measured small alpha changes or the unrestricted projected solution as a qualified repair. Fresh radial validation and Cartesian candidate verification remain on hold.

## Actual results

Every row changes one metal polarizability; all six fitted parameters, geometries and provisional references remain fixed. The opposite metal is reused exactly. No row was selected as a replacement.

| Probe | Weighted residual loss | Passed checks / 25 | W8→W9 paired energy error (kcal/mol) |
|---|---:|---:|---:|
| Baseline | 0.895803724 | 22 | −4.338284 |
| La_minus5 | 1.035742426 | 20 | -4.371455 |
| La_minus1 | 0.903629566 | 22 | -4.344894 |
| La_plus1 | 0.897622610 | 22 | -4.331685 |
| La_plus5 | 0.997230321 | 20 | -4.305389 |
| Dy_minus5 | 0.970072547 | 22 | -4.317806 |
| Dy_minus1 | 0.896591312 | 22 | -4.334204 |
| Dy_plus1 | 0.901807256 | 22 | -4.342357 |
| Dy_plus5 | 0.991443982 | 20 | -4.358582 |

The individual energy limit remains 1 kcal/mol; each rigid-force control requires RMS ≤2 and maximum component ≤5 kcal/mol/Å. No tolerance changed.

## What the response calculation means

Adding two relative-alpha columns increases the retained rank from six to eight. Their finite-difference columns agree between 1% and 5% to within 0.072% in relative norm. The two new directions are numerically resolved under the declared relative singular cutoff 1e−6.

The unconstrained 1% projection reduces weighted residual loss from 0.8958 to 0.1415, but requires La alpha −0.1769 Å³, Dy alpha −4.3616 Å³ and Dy O-pair damping −0.3185. It exceeds the six-parameter bounds, the local trust box and the measured alpha range. Even this inadmissible extrapolation passes only 23/25 predicted checks; Dy/W9 and paired/W9 maximum-force checks still fail. The 5% derivative gives the same conclusion. These are unevaluated linear predictions, not molecular results.

In this native model alpha changes both induced susceptibility and Thole screening length. The test does not isolate either contribution or establish a unique missing physical term.

## Integrity, cost and scope

Job17695 completed0:0 in16s on8CPUs:128 allocated CPU-seconds. It executed exactly40new physical cells/80native programs and reused40 matrix entries from10unique admitted selected cells. Application wall time14.702s; measured scheduler total CPU53.894s. No new quantum calculation, Cartesian finite difference, geometry search or optimizer ran.

Independent raw/input/producer review passes607checks; independent saved projection review passes67checks. Analysis source controls191/191 and native source controls239/239 passed before their respective execution. Raw review, scientific failure and mathematical extrapolation remain separate. Source preparation and review costs are retained in their individual receipts.

The provisional ECP-PBE0-D4 reference uncertainty remains unresolved for a validated exchange prediction. These consumed gas-phase cluster diagnostics do not establish bulk hydration, protein affinity or a bound-minus-aqueous free energy. Donor17621 continues with its original water model; no silent combination occurs.

## Exact evidence

- Native result: `run_v1/execution/RESULT.json`
- Independent native review: `actual_review_v1/REVIEW.json`
- Saved analysis and all predicted controls: `analysis_v1/RESULT.json`
- Independent matrix review: `analysis_actual_review_v1/REVIEW.json`
- Root admission: `run_v1/ROOT_FINAL_ADMISSION.json`

All workspace paths are relative to the matching `workspaces/.../metal_alpha_response_v1` directory.
