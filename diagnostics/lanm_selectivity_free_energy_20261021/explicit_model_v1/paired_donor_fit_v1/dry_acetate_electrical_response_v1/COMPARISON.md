# Dry mono-acetate electrical-response comparison

**The frozen 17621 donor model underestimates the central projected moment and overstates its individual rotation slopes relative to consumed PBE0.** Its Dy−La slope difference has the same sign but smaller magnitude. All 40 model calls and 31 numerical checks pass independent saved review; these checks do not establish physical accuracy.

## Common-origin moments and rotation slopes

Moments are projected along the frozen laboratory tangent. Model values below use field h=.0005 au; the other declared field size is .001 au. Both charged systems have Q = +2: model field energy uses lab 0, and saved quantum moments use μlab=μnative+2Rorigin. No donor u=v assumption enters the energy derivative.

| Central quantity | Model | PBE0 |
|---|---:|---:|
| La μ, ea0 | 0.114374815 | 0.580438037 |
| Dy μ, ea0 | 0.100683834 | 0.635663281 |
| Dy−La μ, ea0 | -0.013690981 | 0.055225244 |

For the original centroid rotation q=(1 Å)θ, dμ/dq=[μ(+q)−μ(−q)]/(2q), in ea0/Å:

| Quantity | q, Å | Model | PBE0 |
|---|---:|---:|---:|
| La | 0.01 | -1.150959515 | -0.718668978 |
| La | 0.005 | -1.150972264 | -0.717913785 |
| Dy | 0.01 | -1.173905782 | -0.789075578 |
| Dy | 0.005 | -1.173915344 | -0.788828538 |
| Dy−La | 0.01 | -0.022946267 | -0.070406600 |
| Dy−La | 0.005 | -0.022943080 | -0.070914752 |

## Numerical qualification and interpretation

Across all ten geometries, halving the field changes μ by at most 7.97×10⁻⁹ ea0 and model χ by 6.38×10⁻⁵ a0³. The unchanged gates are 1×10⁻⁴ ea0 for μ and 0.5 a0³ for individual χ (1 a0³ paired). Model dμ/dq changes by at most 7.97×10⁻⁷ ea0/Å between field steps. All dependent derivative-envelope checks pass. No quantum χ was computed in this panel.

At q=.005 Å and h=.0005 au, model derivative printing bounds are 3.19×10⁻⁶ ea0/Å individually and 6.37×10⁻⁶ paired. Including saved quantum printing/origin precision gives model−PBE0 bounds 3.13×10⁻⁴ individually and 6.26×10⁻⁴ paired. These are rounding bounds, not total numerical or physical error bars. Geometry-step differences remain descriptive; no new tolerance or extrapolation is applied. Native analyze E prints no induction residual, so only the unchanged declared POLAR-EPS and absence of native failure warnings are asserted.

A field-independent vdW adjustment can change mechanical energy/force slopes but cannot change these fixed-geometry μ, χ or dμ/dq values. The discrepancies therefore constrain the electrical part of the frozen model as well as its mechanical behavior. They do not uniquely identify charge transfer, induction, or a particular parameter repair. The centroid rotation also changes metal–oxygen distances; this is not a pure angular test.

## Cost and evidence

**Job 18764: 26 scheduler seconds × 1 CPU = 26 allocated CPU-s.** Application 22.845s; summed native elapsed 18.599s. This comparison adds zero native calls.

[Producer result](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/dry_acetate_electrical_response_v1/run_v1/execution/RESULT.json) (`71d83185`), [independent admission](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/dry_acetate_electrical_response_v1/independent_actual_review_v1/REVIEW.json) (`f21ff157`), [review report](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/dry_acetate_electrical_response_v1/independent_actual_review_v1/REPORT.md) and [independent raw Decimal replay](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/dry_acetate_electrical_response_v1/independent_actual_review_v1/RAW_REPLAY.json) retain both field/geometry steps, all 40 outputs, ten reused zeros, 31 checks and exact source/parameter lineage.
