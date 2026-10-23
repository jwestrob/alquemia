# One constrained local proposal at saved origin 24

Source-only declaration, 22 October 2026. No solve or molecular evaluation is
performed by preparation. The root-owned release must bind this exact source,
inputs and specification before the one numerical solve.

## Fixed information

Use the last original derivative origin
`e63e7189dbecdb8e4112aceb5e4515f69644dbfa574b5fe07644d01ae736a037`,
its exact 24 saved probes and actual normalized spans. This origin passes all
25 original E/F gates. It is not the selected final candidate. Reuse the admitted
673-vector component diagnostic, original domain and literal P8/B8/P9/B9 targets.
No new reference, geometry, parameter evaluation or fitted target is introduced.

The original four residual blocks retain the exact 531-row ordering and half
weights. The component blocks are `[P8,P9]/sqrt(2)` and `[B8,B9]/sqrt(2)` in
kcal/mol with the previously declared 1 kcal/mol engineering scale. Consequently
their squared norms are J_P and J_B. Fixed scaling cancels from the ratios below.

## Frozen proposal problem

Let z be the 12 normalized parameters, d their increment, r_k and J_k the saved
residuals and finite-difference Jacobians at origin 24. One proposal maximizes
`t`, with `0 <= t <= 1`, subject to:

- `||r_P + J_P d||² <= (1-t)||r_P||²` and the same for B;
- `||r_k + J_k d||² <= ||r_k||²` for each original E/F, pair-dipole,
  full-shell-dipole and susceptibility block;
- all ten original energy gates `|e_i + j_i d| <= 1 kcal/mol`;
- all fifteen original body-force gates: category-balanced rigid RMS at most
  2 kcal/mol/Å and every signed body component at most 5 kcal/mol/Å in magnitude;
- `max(-.025,-z_j) <= d_j <= min(.025,1-z_j)`.

There is no new scalar weighted fit loss, tradeoff weight, margin or selection
among starts. The objective is a common fractional reduction of the **squared**
P/B norms; it is not the tangent fractional RMS derivative used in the first LP.
All original block constraints retain complete linearized residual norms.
The constraints are convex quadratic or affine in `(d,t)`. The feasible zero
step is the single initial point. No multistart, trust expansion, fallback solver
or repeated solve is authorized by this declaration.

## Exact individual-gate construction

Read the ten physical energy errors and the fifteen raw body-residual vectors
from each original stored score. Form their Jacobians with the same actual
normalized probe spans. Paired radial work is explicitly retained even though
it is derived and is absent as an independent row from the 531-vector loss.

For each water fixture, raw body residuals are ordered ion net xyz, then water
net xyz and centroid torque/1 Å xyz in the original helper order. The original
category weights are `1/6` on the ion's three components and `1/(12*Nwater)` on
each water component. Thus `RMS²=sum(weight*raw²)`. Reconstruct every original
RMS/max gate and all original 413 residual rows at the origin and probes. At
fixed coordinates these body projections are linear in Cartesian forces, so
finite-differencing the physical components is the same linear model as
projecting the force Jacobians. This is not a prediction of nonlinear native
forces at a new point.

## Solver and prospective numerical admission

Use installed SciPy SLSQP, analytic objective and constraint Jacobians, one
numerical thread, one call, `maxiter=200`, `ftol=1e-12`, and no runtime install.
The returned result is reported with status/message, iteration/evaluation counts,
normalized and physical direction, parameter bounds, all 25 predicted gates,
all six predicted block scores, and constraints at zero, half and full step.
Independently recompute every constraint and box violation. Numerical admission
requires solver success and maximum normalized violation at most `1e-8`; the
raw gate values and any excess over their unchanged physical thresholds remain
explicit. This allowance addresses arithmetic only. No global-optimum or
nonlinear-feasibility claim follows from SLSQP success.

A block norm at most `1e-12` is recorded as a degeneracy and stops this one
proposal without a strict-descent claim or substitute objective. A failed solver
or failed numerical feasibility check is retained; no automatic retry occurs.
A zero solution does not establish failure of the nonlinear model family.

## Interpretation and next boundary

This is a local secant model at one already consumed origin. Even a numerically
admitted nonzero direction is only an unevaluated proposal. The native nonlinear
energy/forces, all original gates, response and positive-definiteness guard must
be checked in a separately released finite evaluation before judging it.
The existing tangent-LP direction, whose full linearized E/F norm grows strongly,
is not submitted for native evaluation. No fresh validation, fit continuation,
new physical term, reference revision or production change is released here.
Record source/input hashes, wall/CPU cost and numerical versions; native and
reference calls remain zero. Root owns any solve and subsequent release.
