# Local linear feasibility of the unchanged twenty-five gates

Declared before calculation. Use only immutable terminal17676 selected data,
its final Jacobian at that exact selected point, raw review17683 and the saved
independent residual review. This tests feasibility of the engineering gates,
not the original least-squares objective and not the physical parameter family.
No native calls, fitting, new outcomes, candidate update or automatic evaluation.

## Fixed local model and necessary linear program

Reconstruct the413×6 Jacobian from twelve saved derivative residuals and their
actual normalized spans. Its point must equal the frozen selected parameters.
Invert the declared residual weights to recover eight raw energy rows and405
raw body-component rows in fifteen groups. Derive both paired radial works as
Dy work minus La work, retaining all ten energy gates and fifteen force gates.
Verify the reconstruction against every saved energy and force component/RMS.

Let delta denote six parameter changes normalized by the original box widths.
Intersect each original parameter bound with the declared trust interval
−0.025≤delta≤0.025. Minimize a common nonnegative threshold multiplier t subject
to every predicted |energy error|≤t×1kcal/mol and every predicted |body component|
≤t×5kcal/mol/Å. This linear program is a necessary relaxation: it omits the
fifteen category-balanced RMS≤2 constraints. A t≤1 result alone is not a pass.

Use SciPy HiGHS dual simplex, one native solver thread, max10000 iterations,
30s solver limit, primal and dual feasibility tolerances1e−9. Report solver
status, optimal objective, primal residuals, dual stationarity/complementarity,
dual objective/gap and every projected gate. Numerical certification tolerances
are1e−7 for normalized LP residuals and dual gap, without relaxing the1/2/5
scientific thresholds. If t lies within1e−7 of1, label the feasibility boundary
numerically unresolved unless an explicit predicted point passes exact gates.

## Conditional norm-constrained solve

If the necessary LP has t≤1+1e−7 but its projected point fails any RMS gate,
perform exactly one SLSQP minimax solve including all fifteen convex norms.
For each global weighted force segment v, RMS/2 = sqrt(20)×||v||₂. The common
t therefore also bounds these fifteen ratios. Start at the LP delta, with t
raised to its largest actual normalized gate ratio. Use analytic Jacobians,
max300 iterations, ftol1e−10 and a30s callback wall limit, one thread. Preserve
failure/status/constraint residuals and the final point. No retries or tuning.

If the LP already supplies a point meeting every original gate, skip SLSQP.
If a numerically certified necessary LP optimum exceeds1+1e−7, skip SLSQP and
report infeasibility only of this local trust-box linearization. Solver failure
or an unsuccessful norm solve is not a certificate of physical-family failure.

Every solution remains a prediction from the saved local model. Preserve all
twenty-five gates and disclose that the Jacobian uses finite normalized halfsteps
0.005, the trust interval is finite, and the residuals are correlated engineering
controls. Measure one-thread wall/CPU costs and preserve all parent bytes.
