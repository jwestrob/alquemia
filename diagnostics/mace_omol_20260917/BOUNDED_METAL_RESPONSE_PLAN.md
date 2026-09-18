# Bounded metal response — next contained development test

Declared2026-09-18 after completing METAL_RESPONSE_PLAN.md, under Jacob's
active autonomous discriminator goal. This is a NEW descriptor protocol;
`DFT_anchored_MACE_GB_full_metal_translation_v1` and its failed alpha eligibility
remain unchanged. All four original representations are consumed development.
No new benchmark label, calibration threshold, water state or protonation.

## Question and cost

The full3D model passes all four native GGR energy/gradient checks and gives
response corrections1.189859/1.365836kcal for extended/connected cores; partition
change0.175977kcal. Alpha's free quadratic optima lie0.244–0.349A from their
centers, outside the predeclared0.20A domain. The next question is whether a
fixed, limited amount of geometry response helps discrimination, even when
an unconstrained stationary minimum is unavailable in that domain.

Apply the SAME rule to all four representations and both metals. Reuse every
existing center, cheap grid and gradient. The GGR solutions remain exactly the
already executed interior points and should reuse their four nativeDFT and
eight short records after exact physical/method matching. Only four new alpha
analyticDFT endpoints and eight matching core/full short evaluations are needed.
No new MACE-core/GB grids, numericalDFT derivatives, trajectory or quantum
optimization. This remains two high-level endpoint positions per state.

Protocol ID: `DFT_anchored_MACE_GB_bounded_metal_response_v1`.

## Fixed finite-response rule

Use the same half-step3x3 matrices K=K_C+K_J and center gradient g=g_DFT+g_J.
The radius remains r=0.20A, frozen in the original experiment before outputs.
For every state solve the positive-curvature quadratic problem

    minimize q(u)=g^T u + 0.5 u^T K u, subject to ||u||_2 <= r.

If the unconstrained solution lies inside the sphere, preserve its recorded
coordinates exactly. Otherwise solve

    (K + lambda I)u = -g, lambda>=0, ||u||=r.

Use a deterministic bracket/doubling and bisection of lambda (at most200 tiny
3x3 solves; norm residual<=1e-10A and equation relative residual<=1e-10).
This is the actual quadratic minimum on a sphere, not componentwise clipping
or eigenvalue clamping. Preserve all eigenvalues and reject nonpositive coarse
or fine K. Predict the energy using g^T u+0.5u^T K u; the unconstrained
-0.5g^T solve(K,g) expression is invalid for an active boundary.

Keep the original numerical convergence conditions: C,J,combined curvature
refinement changes at the new point<=0.02kcal, and coarse/fine displacement
change<=0.01A. Require identical source/caps/protonation/waters/donor inventory.
No changed donor cutoff or selective coordinate rescue. Failed states remain
unavailable. Radius is a protocol parameter, not a fitted biological label.

The scored descriptor is the ACTUAL evaluated DFT+J change at this fixed point,
not an assertion that the unconstrained physical system has equilibrated:

    B_M = E_DFT,M(u_M)-E_DFT,M(0) + J_M(u_M)-J_M(0)
    R_bounded = R_center + B_Ca-B_La.

No second full solvent energy, entropy, aquo redefinition, old decision bands,
or combination with the rejected conductor correction. It omits long-range
scaffold electrostatics just as the parent conditional mechanics model does.

## Independent native checks and usefulness

Freeze all proposed points before newDFT outputs. Preserve the original energy
check: prediction error<=max(0.05kcal,25% of the nonlinear DFT energy change),
and actual energy lower than the center. Interior points retain the original
full-gradient criterion<=max(2kcal/mol/A,25% of initial gradient norm).

At an active sphere boundary, a nonzero radial gradient is expected. Require
its radial sign consistent with the constraint (u·g_actual<=0) and its tangent
component norm<=max(2kcal/mol/A,25% of initial gradient norm). Report the entire
gradient and radial component; do not call this unconstrained stationarity.
This is energy/gradient validation at a fixed displacement, not a nativeDFT
Hessian or global-basin proof. The same GGR partition requirement<=2kcal applies.

Report actual exploratory contrasts even if a physical check fails, with final
qualified fields separate. Primary usefulness remains both alpha structures
above GGR under the shared extended-core policy, without refitting. Report
margin changes separately from any successful direction changes. Partial
improvement does not establish broad discrimination. Structural replicates
remain grouped and no previously consumed case becomes a blind test.

## Status

Plan prepared; no new bounded-protocol implementation or calculation has run.
The parent unconstrained pilot is complete and will retain its exact results.
Proceed through the existing runners and pinned environments, keeping the
baseline/default and other agents' work intact. Record all reuse and actual
costs. No automatic promotion or production rescore.
