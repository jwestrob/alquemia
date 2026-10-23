# One saved ceiling-tradeoff diagnostic

This is a separate diagnostic. The original full conic result remains AlmostSolved, its full-fraction epigraph failure remains recorded, and `conic_joint_ceiling_removal_v1` remains unrun under its unchanged Solved gate. No fit or physical threshold changes.

## Baseline certificate, before any new solve

Use the independently reviewed saved origin-24 conic matrix, physical model and full returned twelve normalized increments. Recompute the original CORE scalar certificate at fraction one. Define rho as the maximum of the two recomputed P/B squared-error ratios. Freeze auxiliary q = rho + 1e-8 * max(1, abs(rho)); this is epigraph headroom only, never a physical parameter or feasibility tolerance. Require all 851 original residual margins strictly positive, all 25 physical gates, all four original score ceilings and exact parameter/trust bounds to pass unchanged. Preserve the unmodified raw q and its failed certificate. No solver is needed for this certificate. Do not alter a parameter to obtain a pass.

## One prospective solve

Reuse exact saved affine arrays and delete only the four development ceiling cones (535 rows). Retain both component epigraph cones, all 25 physical gates, every parameter bound and the normalized +/-0.025 trust region at the original origin. Keep the original minimax objective t >= 0, q=t^2. The resulting matrix has 1282 rows, 13 columns, 17 SOCs and 856 nonnegative sides. No native evaluator or new reference is imported.

Use every original Clarabel 0.11.1 setting, including 200 iterations, 60 seconds, one thread, qdldl and 1e-10 primary tolerances. Reserve one persistent exclusive attempt before constructing the solver; one solve only, no retry or fallback. Root release and exact source review are required before that attempt. Baseline preparation is authorized separately and calls no solver.

Preserve raw solver status, x/s/z, all residuals, dual gap and cone margins before optional reporting. Status is not a witness-feasibility test. For each original fraction 1, 1/2, 1/4, 1/8, interpolate q as CORE does; independently require all 847 retained scalar margins strictly positive, all 25 physical gates and bounds. This is ordinary-double numerical feasibility, not formal interval proof. Always report the four removed scores and whether their original ceilings would pass. No physical candidate is selected or admitted.

For an independently feasible witness compare its actual max(P/B squared-error ratios), not its auxiliary q, with the baseline's actual rho. Report signed gain and the inherited 1e-10*max(1,abs(rho)) numerical comparison flag. The baseline cushion must not create apparent gain. A gain beyond the supported full-problem objective uncertainty supports a local preservation tradeoff; failure to find a gain cannot establish infeasibility or missing physics without adequate dual/optimality evidence. AlmostSolved stays AlmostSolved. This test says nothing about the nonlinear model outside the fixed local span.

All authorities are immutable pins, including the full actual arithmetic review and stored sparse arrays. Output and attempt marker paths are fixed in INPUTS. Resource envelope: one numerical thread; zero native/reference/heldout calls. No change to current campaign, parameters, scientific tolerances, original results or the old gated diagnostic.
