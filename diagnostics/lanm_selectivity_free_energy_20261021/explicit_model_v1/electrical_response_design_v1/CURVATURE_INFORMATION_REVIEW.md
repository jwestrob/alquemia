# What a finite-field experiment would newly measure

Source/design review after the admitted zero-field trimer and 18413 sensitivity results. No implementation, field amplitude, molecular execution, submission or engine change is made here. The nonzero pair-response prototype still rejects external fields.

## New information: curvature, not the known slope

For one fixed laboratory direction e and scalar field f in consistent units,

`E(f) = E0 − mu0 f − alpha f²/2 − beta f³/6 − gamma f⁴/24 + ...`.

The saved seven-term E3 and mu3 already give the value and first derivative of that same cycle at zero field: `dE3/df = −e·mu3`. Merely confirming the slope would repeat existing information. The new observable is `alpha = −d²E/df²`, or a declared balanced combination of these curvatures, with a check for resolved field dependence.

For the baseline linear induced-dipole model at fixed geometry, `A = alpha_site_inverse − T_mutual` and `u(f) = A^-1(E0 + B f)`, where B applies the uniform field to every polarizable site. Its molecular susceptibility is `B^T A^-1 B`, with the native unit factors restored. This is independent of the permanent driving field E0. Thus a direct-only kernel change can alter zero-field energy/dipole without changing uniform-field susceptibility; a mutual-kernel change can change susceptibility. This is an algebraic distinction, not permission to bypass the prototype's field guard.

Native `tinker_source/source/induce.f:220–225` adds the same external field to both driving channels after their permanent fields are computed. The existing [CPU field source review](REPORT.md) records the correct field-energy route, unit conversion and the scoped water-group equality of the two channels. The dual stationary functional and the distinction between its direct and mutual kernels remain as documented in `collective_mechanics_options_v1/DUAL_FIELD_IMPLEMENTATION_MAP.md`.

## Smallest useful conditional comparison

A four-trimer near/far × La/Dy field panel would measure new total susceptibilities, but would not isolate the known collective defect from pair or individual-water response. At the same 16-reference-point cost, a more directly matched conditional question is:

**Does adding the second water to one exact consumed M–water dimer change directional susceptibility by the same amount in baseline mechanics and provisional PBE0?**

Use one predeclared geometry-only context, for example the existing W9 near pair (waters 1 and 5), with its exact M–water-1 fragment, for both La and Dy. Preserve parent coordinates, full-parent quantum ghost basis, state, and all existing model parameters. Apply one shared laboratory direction to each trimer and its dimer; define it from geometry before outputs, without separately rotating fragments. For every system, retain `+h, −h, +h/2, −h/2` and reuse its admitted zero-field energy/dipole.

The observable is `Delta alpha_env(M) = alpha(MW1W5) − alpha(MW1)` and its Dy−La difference. Four systems × four signed fields = **16 new reference SPs and 16 baseline model analyze calls**, before any separately required native qualification controls. Both terms have charge +3, so their energy/dipole difference is charge neutral in the unchanged parent frame. The field curvatures of each fixed-charge system are themselves origin invariant. This test needs neither new geometry nor a new quantum method.

This increment includes the added water's intrinsic susceptibility, its pair response, and the response of the pre-existing fragment. It is **not alpha3**, not a unique collective-metal susceptibility, and not unique charge-transfer or polarization attribution. Comparing model-minus-QM increments tests an environment-dependent response error, but cannot assign its microscopic source. If root requires an explicitly irreducible three-body curvature, the smaller conditional panel does not answer that stricter question.

The exact cost of that stricter question is material: one seven-term trimer cycle for both metals at four fields needs 7×2×4 = **56 reference SPs**. For W9 near(1,5)−far(1,3), M, W1 and MW1 cancel within each metal's identical parent basis. Eight remaining fragments per metal give **64 SPs**. Neutral-water fragment fields do not cancel exactly across La/Dy because their ghost metal bases differ. Rotational replication would introduce an unqualified numerical assumption. These counts were cross-checked with Singer; they are scope distinctions, not additional proposed runs.

## Required curvature and field-step controls

Use both energy and moment readbacks from the same native calls:

- `alpha_E(h) = −[E(h)+E(−h)−2E0]/h²`.
- `alpha_mu(h) = [mu_e(h)−mu_e(−h)]/(2h)`.
- Repeat both estimators at h/2. For smooth response, their leading quartic corrections are respectively `gamma*h²/12` and `gamma*h²/6`. The combination `[4*alpha_E(h/2)−alpha_E(h)]/3` removes that leading energy error, but is a diagnostic estimate, not a substitute for checking numerical noise and the two raw scales.
- Preserve the odd energy part and even dipole shift against the already known mu0. They check field sign/units and can expose higher response; they are not the primary new scientific outcome.

Freeze h, direction and acceptance/error budgets before field outputs. Require energy-curvature/moment-slope consistency at both steps, exact method/geometry/state and site denominators, finite native residuals, and unchanged electronic-state admission. Retain every signed result. A half-step discrepancy alone is not evidence of physical nonlinearity: printing, SCF branch changes and incomplete induction can mimic it. No automatic step extension is implied by a failed precision check.

The existing model energy rounding bound is 5e−9 kcal/mol per value. A single-system energy-curvature numerator therefore has bound 2e−8 kcal/mol; balanced coefficients add their absolute rounding bounds. Dividing by (h/2)² amplifies the numerator bound fourfold relative to h. Six-decimal saved induced dipoles must supply model moment slopes; three-decimal printed total moments are only a consistency check. Their fixed permanent and COM-origin contributions cancel between opposite fields of the same unchanged fragment. Solver residual thresholds alone are not proven observable-error bounds. Reuse the explicit native/QM unit conversion in the earlier source review.

At fixed geometry the baseline model is linear in field and has an exactly quadratic field energy before numerical error. Resolved classical nonlinearity would therefore be an implementation/admission concern. Resolved smooth QM nonlinearity would be genuinely new evidence beyond its zero-field moments, while still requiring a stable electronic branch and sufficient precision.

## Interpretation boundary

- Agreement in the conditional curvature despite a wrong zero-field dipole increment would favor a defect in zero-field driving/offset terms for this measured mode over adjusting its susceptibility. It would not prove that the whole response matrix is correct.
- A significant curvature discrepancy would show that the baseline response in this environment and direction is wrong in addition to the known offset. It would support investigating the response operator or its environmental dependence, without uniquely identifying charge transfer, penetration, exchange, or a particular parameter.
- A uniform field probes only spatially uniform modes. The metal's internal field is inhomogeneous; one direction does not reconstruct the tensor, local response kernel, or a bulk/donor response. Even a complete tensor would not uniquely determine those internal modes.

The completed E3/mu3 and 18413 evidence remains unchanged. No field result would erase the original full coordination-energy failure or qualify provisional PBE0 as an exact reference. This review supplies a distinct information test, not another fit or a physical repair claim.
