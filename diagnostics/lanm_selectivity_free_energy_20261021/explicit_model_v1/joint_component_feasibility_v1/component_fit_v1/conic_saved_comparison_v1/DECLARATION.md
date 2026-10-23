# One saved-origin solver comparison

Question: can a solver designed for convex cone constraints solve the exact
already consumed origin24 linearized proposal with a useful primal/dual residual
and gap, without changing any physical gate or collecting another molecular cell?

This separate engineering diagnostic must not modify, replace or extend the
running eight-step campaign. Execute only after that campaign is terminal.
Use only its immutable origin24 and exact first24 saved derivative vectors,
original CORE linearization and fixed ceilings/denominators/domain. Preserve the
original SLSQP result; do not rerun SLSQP or select an origin from later outcomes.
Use the independently reviewed conic mapping(b75d232e) in exact correspondence
with CORE. All residual gates, normalization and trust bounds remain unchanged.

Prepare Clarabel0.11.1 in an isolated environment, retaining current numpy/scipy
versions as read-only system packages. Never install into or change the active
fit environment. Pin installed versions, wheel metadata and exact settings.
Primary interfaces/settings: https://clarabel.org/stable/python/getting_started_py/
and https://clarabel.org/stable/api_settings/ ; package identity from
https://pypi.org/pypi/clarabel/0.11.1/json .

One actual solver call, maximum200 iterations/60seconds, one thread, qdldl,
absolute/relative gap and primal/dual feasibility tolerances1e-10; infeasibility
tolerances1e-10. Freeze all remaining effective settings before that solve.
Retain returned status including reduced-accuracy status, raw primal/dual vectors,
gap and residuals. Compare q=t^2 and the old four fractions through the unchanged
independent scalar certificate; interpolate in q, not t. No new fraction search,
native candidate evaluation, physical qualification or production change.

Source controls may check cone/CORE constraint correspondence on tiny artificial
problems and fixed rational probes, without an actual scientific solver call.
Source review and metadata pins precede actual execution. Report mathematical
feasibility separately from optimizer convergence and strict improvement.
Floating solver dual residuals are numerical evidence, not a rigorous proof of
global physical infeasibility or nonlinear optimality.

Maximum new native/reference/held-out calls:0. Maximum actual proposal solves:1.
No automatic subsequent fit is released by this diagnostic.
