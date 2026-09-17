# Numerical convention fixed before cheap grid execution

The immutable PLAN.md specifies coarse/half-step matrices but does not choose
which supplies the response prediction. Freeze this choice now, before the new
MACE/short grids run or any new DFT energies are inspected: use the HALF-STEP
matrix for every response prediction, conditional on the prescribed coarse/fine
convergence gate. Both matrices remain reported; no per-case choice is allowed.
DFT curvature/anchored checks compare this chosen fine matrix against the actual
coarse displacement energies. This also tests finite-amplitude approximation.

Use dimensionless coordinates z=(s/h_s,theta/h_theta), where h_s=0.02Angstrom and
h_theta=pi/180radian. For a grid spacing a (1 or0.5), compute

    A_ii = [E(+a e_i)+E(-a e_i)-2E(0)] / a^2
    A_12 = [E(+a,+a)-E(+a,-a)-E(-a,+a)+E(-a,-a)] / (4a^2).

Here A=diag(h) K diag(h) has energy units. The physical gradient transforms as
b=diag(h)g, so solve A z*=-b and deltaE=-0.5 b^T solve(A,b). Do not mix Angstrom
and radians in an unscaled eigenproblem or a trust test. Eigenvalues refer to
this fixed dimensionless measure; positive definiteness is coordinate invariant
under this nonsingular scaling. Do not interpret these as mass-weighted modes.

For the convergence gate, find the exact largest absolute value of
0.5 z^T(A_half-A_full)z on [-1,1]^2 using all corners, stationary points on each
edge, and zero (all interior stationary points of a homogeneous quadratic have
zero value). Preserve the mixed term and signs. Missing/failed matrices are
unavailable, never replaced by zero or the coarser matrix.

This note adds a deterministic numerical convention; physical inputs, task
counts, energy expression, trust region, tolerances and hypotheses are unchanged.

At an eligible, fixed predicted minimum, retain the quadratic relaxation as a
prediction and compare it against the single actual DFT+J energy change. If the
validation gate passes, use that ACTUAL evaluated change as the primary response
descriptor, with the quadratic estimate reported separately. This uses the
second high-level endpoint per state and avoids discarding its information.
The position is fixed before that output; there is no follow-up optimization or
choice of whichever energy gives a desired class. Numerical positivity must
hold for both fine and coarse combined matrices; report all eigenvalues without
clipping. These choices are frozen before cheap grid execution or result inspection.
