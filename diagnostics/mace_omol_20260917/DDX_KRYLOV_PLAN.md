# Solve the same native PCM equations with a maintained Krylov method

The unchanged native Jacobi/DIIS solves remain slow across the finer real
protein representations. One coarse paired result also fails reciprocity.
Numerical iteration and representation errors must be distinguished. This
engineering experiment changes the linear solver, not the dielectric model.
Keep all ongoing jobs and their frozen implementations untouched.

## Native operator access and isolated build

Add a minimal, separately built adapter exposing the existing ddX0.9.0
`rinfx`, `repsx`, `lx`, `prec_repsx`, `ldm1x` and `hnorm` routines, plus the
source vector already assembled by native `pcm_setup`. Do not replace those
routines, discretization, physical parameters or source construction.
Native Jacobi uses off-diagonal operators: a full operator call must set
`constants%dodiag=.true.` and restore the original flag afterward. The two
native equations are `R_epsilon * phieps = R_infinity * phi_native` and
`L * x = phieps`; `phi_native` already carries the native negative source sign.
Do not introduce another sign reversal.

Build a new isolated wheel/environment from the pinned source archive and
existing pinned dependency wheels. Preserve original and patched source hashes,
the exact patch, build logs, imports, compiler/library versions and costs.
No working environment or current executable may be replaced. Build/import
alone is not a scientific qualification.

## First real qualification: two existing coarse endpoints only

Use GGR2FW0 coarse Ca and La, both with actual completed native solutions from
the recovered source group. Preserve the exact source/cavity,6/194basis/grid,
FMM12/12,eps78.3,eta.1,shift0,64threads and matrix-free storage. Two new forward
endpoint attempts, at most four linear-system solves, no automatic retry.
Use SciPy1.17.1 GMRES, restart40,maxiter30restart cycles,rtol1e-10,atol0,
native block preconditioners, and zero initial guesses. Solver iteration
settings are numerical controls, not a project compute/time budget.

Transform coordinates using weights `1/sqrt((l+1)*n_spheres)` so the Euclidean
norm is the native H^(-1/2) norm. Freeze these checks before execution:

- Transformed norm versus native hnorm: relative error<=1e-12.
- Independently reevaluated full-operator residual in each system:
  relative native hnorm<=1e-9. Record both GMRES's return status and residual.
- Saved native solutions: evaluate the composed equation
  `R_epsilon*(L*x)-R_infinity*phi_native`, without a new native iterative solve;
  report its relative residual with the same1e-9 screen.
- Each endpoint and Ca-minus-La energy agrees with the actual native solutions
  within1e-6kcal. Energy remains0.5*dot(psi,x), converted once.
- Exact cavity-point identity, source phi/integral checks and passive-energy
  checks remain unchanged. Preserve all failures, zero/missing values correctly.

Report setup, operator/preconditioner counts, iteration history, both true
residuals, energy differences, wall/CPU/memory and allocated costs. This cannot
repair a coarse discretization reciprocity failure; retain that failure. No
DFT/MACE/gradient/biological score, new physical state or calibrated threshold.
No wider Krylov rescore is included in this first qualification. An extension
requires a new recorded scope, executed autonomously under the active goal.
