# Mechanics design recommendation — source/evidence only

**Recommend testing independent radial drive shape before changing angular response or adding a generic collective energy term.** This is a discriminating hypothesis, not an identified missing physical mechanism. No source, coefficient, fit or numerical calculation was changed for this recommendation.

## What the current family already contains

The [functional](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/finite_screen_v1/interpretation_v1/SPEC.md) and [native_helper](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/native_v1/LOCAL_RESPONSE.f) implement, for each metal,

`h_i = 3 f(r_i) phi(r_i) r_i [a + b eta_i]`,
`H_MM = [exp(g(s)(c+d eta))-1] I/alpha_M`.

The eight response amplitudes independently control pair/collective drive and pair/collective isotropic metal stiffness. Four additional parameters are native buffered-14-7 metal–O R/epsilon. Thus coordination-dependent response and distinct driving/stiffness channels are already present. Water–water interactions, intrinsic polarizabilities, permanent multipoles and native mutual T remain fixed.

| Freedom | Current restriction or untested alternative |
|---|---|
| Radial value versus slope | f has fixed2 Å⁻¹ shape and2.50 Å normalization; phi uses frozen native damping. a/b/c/d change amplitudes, not an independent radial shape. R/epsilon offer only the native buffered-14-7 shape family. |
| Coordination | eta_i, eta and g are fixed rational functions of radial overlap sums; one collective basis per drive/stiffness. There is no separately adjustable coordination-shape function. |
| Angular | The baseline tensor/multipoles already contain angular physics. The added field is radial on O and the added metal stiffness is isotropic; there is no independent angular descriptor, anisotropic correction tensor or movable permanent-charge distribution. |
| Electrical | h enters both dual RHSs; H acts only on the metal. At fixed geometry, changing drive changes moments/energy but not chi. Changing H changes chi. There are no added water-site stiffness or metal–water off-diagonal response blocks. |

These are structural restrictions, not proof that each requires an extension. The [domain](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/PARAMETER_DOMAIN.json) bounds are development choices; u/v box coordinates must not be confused with native c/d.

## Evidence selecting the first hypothesis

The [actual_two_point](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/two_direction_native_validation_v1/run_v1/execution/RESULT.json) closely reproduces the affine contrast: B-half passes25/25 with P/B RMS12.040/10.192 kcal/mol; D-half reaches6.361/5.193 but reproduces13 failures, dominated by single-water radial mechanics. Removing box faces gave little extra predicted improvement, while preserving three electrical norms cost only a modest amount in the [four_affine_cases](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/final_local_attribution_v1/run_v1/execution/affine_cases/RESULT.json). This points first to the link between interaction-energy values and radial derivatives in the present basis. It does not show that isotropic response is globally impossible.

A purely collective addition that vanishes for an isolated pair cannot directly change P, which is a literal sum of dimers. A field-independent short-range energy patch cannot directly change the W9 dipole sign. Either alone would leave part of the stated question unanswered.

Both central W8/W9 shells in the [prepared_geometry](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/evaluator_v1/preparation_v4/MANIFEST.json) use the same2.50 Å metal–O radius by construction, up to printed-coordinate precision. Changing a decay range while renormalizing f(2.50 Å)=1 leaves their central h/H and energy/moment values unchanged at fixed amplitudes. It supplies derivative freedom, not a new direct central-dipole direction. This distinction is essential to a meaningful test.

## Smallest concrete radial hypothesis

One additional amplitude e_M per metal (two new coefficients total), with a fixed smooth basis and no new fitted length scale:

`psi(r) = (r²-r0²)/r0²`, `r0=2.50 Å`,

`delta h_i = 3 e_M f(r_i) phi(r_i) psi(r_i) r_i`.

Keep f in the coordination sums, all H blocks, T, intrinsic alpha, permanent multipoles and vdW conventions unchanged. Add delta h to both RHSs in the same stationary functional. This makes a field-value/slope basis: at r0 the added field vanishes but its radial derivative does not. Existing a/b and R/epsilon can then move shell energies/moments while the added radial basis compensates some derivative costs. It is **not** a claim that varying e alone fixes P/B or the central W9 sign: those values are unchanged at the ideal central shells when every old parameter is fixed.

The analytic correction uses the existing dual form `−K/2 (u+v)·d(delta h)` plus native force changes through the newly solved u/v; it must include the psi derivative and metal reaction. No extra H energy is added. Contact remains finite and the polynomial-times-exponential field decays at separation. A future domain and numerical implementation still need declaration; neither is supplied or released here.

The smallest discriminating study would compare the current12-column local model with these two declared response-shape columns, retaining all25 physical gates and reporting signed electrical targets. It must ask whether the added columns allow the *old* directions to improve both P/B and W9mu, rather than demand an effect from e alone or immediately run another multistep fit. A negative or rank-deficient local result would reject this radial hypothesis at that origin, not the Hamiltonian class globally.

## Limits and alternative hypothesis

This drive-only extension leaves chi unchanged at fixed old parameters; correcting absolute susceptibility still requires existing c/d changes or a separately justified stiffness-shape extension. It adds no independent angular direction at a fixed-radius shell. If mechanics improves but the W9 dipole sign remains locally unreachable, the next distinct hypothesis is an angularly resolved collective drive or anisotropic response correction that distinguishes ligand arrangements at equal radial coordination—not another amplitude multiplying the same eta. That alternative should be tested separately, since it cannot directly fix P when constructed to vanish in dimers.

Risks are correlation with a/R/epsilon, fitting the narrow consumed radial window, adverse shell force changes (psi is zero there but its derivative is not), and nontransferable response. The present provisional common-basis references and sparse angular sampling cannot identify charge transfer versus polarization versus short-range exchange. No new term, acceptable error threshold, model promotion or sampling readiness follows from this note.

## Exact source pins

- [functional](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/finite_screen_v1/interpretation_v1/SPEC.md): `34cc2c4e03ce667ce033cb91662f86a5603392d7b740079eb9b46473b84e2a65`
- [native_helper](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/native_v1/LOCAL_RESPONSE.f): `16a85d3a55f0e15a3d76469223a6dee15845861de274955a3396e517f206c13a`
- [domain](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/PARAMETER_DOMAIN.json): `a9e53eb59749b45907343e72fe85290e13402c9a20326240861a1028b927799c`
- [prepared_geometry](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/evaluator_v1/preparation_v4/MANIFEST.json): `8f29dde62c7462cfd73960fa237e33f563d51d7c80e3deb1e1c957495042a623`
- [actual_two_point](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/two_direction_native_validation_v1/run_v1/execution/RESULT.json): `1604c8ff0cbaaae2eea3dd7db1f42428826166e39ccc8384e4f9af5f8cee97bb`
- [four_affine_cases](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/final_local_attribution_v1/run_v1/execution/affine_cases/RESULT.json): `346a018d4beabd7e6f8255e55dfbee98dd495935d447d934978808eb0c13b48d`
