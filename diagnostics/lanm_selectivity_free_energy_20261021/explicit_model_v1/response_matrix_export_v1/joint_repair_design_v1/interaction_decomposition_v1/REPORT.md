# Saved pair/collective energy cancellation

The selected candidate's improved W9-minus-W8 coordination contrast contains substantial cancellation between inaccurate pair and collective contributions. Its common-basis coordination error is **−0.3131 kcal/mol = +3.1483 pair −3.4614 collective**. The pair error increased from baseline; the collective error decreased. The individual W8 and W9 Dy-minus-La interaction errors both worsen in this diagnostic. Thus the 25/25 fitted E/F result does not establish accurate pair or collective energies.

## Definition and available evidence

At each of the four unchanged central W8/W9 × La/Dy geometries, use the literal saved cycles:

```
C = E(M+shell) − E(M) − E(shell)
P = Σ_i [E(M+Wi) − E(M) − E(Wi)]
B = C − P
```

C is the full metal–shell interaction; P sums independent metal–water interactions; B is the remainder containing all nonpair orders. B is not a unique three-body or polarization term. The selected46 artifact, baseline/static evaluator terms, all34 pair coefficient dictionaries and admitted common-parent-basis reference terms are pinned in [INPUTS.json](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/interaction_decomposition_v1/INPUTS.json). The consumed eight-trimer result supplies the same-candidate join and motivation.

**Absolute per-parent electronic-plus-D4 targets and their total-energy errors remain unavailable.** The exact real-shell D4 authority is absent; electronic-only C/P/B values are retained in the result without treating them as total-energy error targets. For Dy-minus-La at each shell size, the identical literal water geometries permit cancellation of the unknown own-shell/water D4 terms. Each parent's electronic ghost-basis terms are retained explicitly. Saved full-cluster, dimer and ion D4 terms are then added once. The inherited conversion is 627.509474 kcal/mol per hartree.

## Paired reference components

All entries are kcal/mol. W8/W9 denote Dy-minus-La; coordination denotes W9 minus W8.

| Quantity | Electronic C | D4 C | Electronic P | D4 P | Electronic B | D4 B |
|---|---:|---:|---:|---:|---:|---:|
| W8 | −59.7700 | +0.6748 | −75.1868 | +1.4678 | +15.4168 | −0.7929 |
| W9 | −60.0862 | +0.7815 | −84.6160 | +1.6513 | +24.5298 | −0.8698 |
| Coordination | −0.3161 | +0.1066 | −9.4291 | +0.1835 | +9.1130 | −0.0768 |

## Errors against the available electronic-plus-D4 paired targets

Errors are model minus reference, kcal/mol; each row satisfies C error = P error + B error.

| Quantity | Model | C error | P error | B error |
|---|---|---:|---:|---:|
| W8 | Baseline | +2.7602 | +12.3586 | −9.5984 |
| W8 | Selected | **+4.5473** | **+24.9402** | **−20.3929** |
| W9 | Baseline | −1.9028 | +13.9342 | −15.8370 |
| W9 | Selected | **+4.2343** | **+28.0885** | **−23.8542** |
| Coordination | Baseline | −4.6630 | +1.5756 | −6.2386 |
| Coordination | Selected | **−0.3131** | **+3.1483** | **−3.4614** |

At both shell sizes, the pair contrast is insufficiently Dy-favoring and the collective contrast compensates in the opposite direction. Both component errors grow after fitting. The selected W8 and W9 C errors also become similar (+4.5473 and +4.2343), so subtracting shell sizes removes most of their common error.

The signed coordination C error shifts by +4.3499 kcal/mol toward zero: the B error shifts by +2.7772 and the P error by +1.5727 kcal/mol. Both signed shifts contribute to the smaller total error, although the absolute P error worsens while the absolute B error improves. Its small final residual therefore does not demonstrate accurate components. The original fitted gates retain their original energy/reference definitions; this common-basis decomposition does not replace those targets.

## All four model changes

Selected minus baseline, kcal/mol. These are model changes, not per-parent reference errors.

| Parent | ΔC | ΔP | ΔB |
|---|---:|---:|---:|
| W8 La | −12.3321 | −13.6633 | +1.3312 |
| W8 Dy | −10.5449 | −1.0817 | −9.4632 |
| W9 La | −10.8136 | −15.3713 | +4.5576 |
| W9 Dy | −4.6765 | −1.2169 | −3.4596 |

La gains binding mainly through the pair sum, partly opposed by a more positive collective remainder. Dy's pair sum changes little; its binding gain mainly comes from reducing the collective remainder. An improved difference of full-shell B values can coexist with worsened sampled E3 energies because B includes all nonpair orders and all water combinations. This calculation does not identify which higher-order terms compensate or establish a unique physical mechanism.

## Limits, provenance and cost

The PBE0/ECP common-basis reference remains provisional. The admitted native records retain their original stopping/state limitations: all68 new and all12 reused terms fail the printed density/DIIS criteria; native Energy Check/Gradient Check endings and other flags remain recorded. These flags are neither erased by cancellation nor converted into numerical error bars. No physical tolerance, fit, fresh-data read or promotion is added.

This diagnostic originally bound the root-admitted selected energy/moment subset for consumed regression; that input scope and its saved result remain unchanged. The full historical fit audit has since completed through the [preserved-phase continuation review](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/fit_v1/independent_actual_review_v1/continuation_v2/actual_v1/REVIEW.json) (`8658d5af`) and [root final admission](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/fit_v1/run_v2/ROOT_FINAL_REVIEW.json) (`1fe38e62`). The audit status is `passed_using_preserved_completed_phases_and_narrow_continuation`; original failed reviews remain preserved. Literal geometry/term joins, the N-fold ion subtraction and the previous common-basis/baseline aggregate were replayed from saved data; arithmetic residuals are below 1.2×10⁻¹³ kcal/mol. No raw native re-admission suite was repeated.

[Result](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/interaction_decomposition_v1/actual_v1/RESULT.json) SHA256 `a9ea13bf6807451de1cb748085df2a8cbfe3686cfadfaa7782993699a277bc57`.
[Declaration](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/interaction_decomposition_v1/DECLARATION.md) and [calculation](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/interaction_decomposition_v1/CALCULATE.py) were frozen before the saved arithmetic. Actual calculation cost: **0.1930 s wall, 0.1613 CPU-s, one thread, zero native calls and zero scheduler allocation**. All four parents and both paired contrasts are retained.
