# Sparse numerical stability check: next implementation

## Recommendation and scope

Use **standard CHOLMOD sparse Cholesky on the shifted real-space matrix**, with the existing positive-semidefinite PME contribution and an upper bound on its norm. This can succeed when block Gershgorin fails because it preserves the signs and cancellation between all real-space matrix entries. First qualify a standalone numerical backend on small matrices and saved vacuum matrices; then separately declare one bounded static-box feasibility run. No periodic engine, trajectory, candidate composition or box factorization is released here.

The inherited guard uses ordinary-double dense eigenvalues. The proposed replacement should likewise be a **qualified numerical sufficient check**, preserving the relative spectral threshold and symmetry/error reporting. A formal interval proof is a different, optional level of assurance; it is not a new mandatory scientific gate. The optional exact-rational tiny control in this directory illustrates that distinction rather than prescribing a bespoke production verifier.

The failed cheap bounds remain closed. Their negative values do not establish that the real matrix, or the complete operator, is indefinite. Sources, library inventory and the 13 tiny artificial controls are pinned in `READY.json`; no box matrix was assembled or factored.

## 1. Exact operator and unchanged threshold

On the actual active positive-alpha induced coordinates, write

\[
 \bar A=B_0+E H_m E^T+R,
 \quad B_0=I-S(T_{real}+T_{self})S,\quad S=\mathrm{diag}(\sqrt\alpha).
\]

Here E selects the metal's three coordinates and H_m is its **normalized** diagonal stiffness. For the present isotropic metal correction, H_m=beta I3 with beta=exp(sigma)-1. The drive h changes the RHS, not this stability operator. All 11,958 sites in these two boxes have positive alpha, giving n=35,874; a general interface must exclude native zero-alpha/nonpolar unknowns while preserving their permanent fields.

The PME source proof gives R=Cbar^T F* Q F Cbar positive semidefinite only with the specified current orthorhombic/conducting initialization, nonnegative conjugate-even multipliers, and matching scatter/gather. Do not manufacture evenness by symmetrizing native multipliers. Preserve the original 0.4/A and 9 A engineering split; bind the eventual actual grid/spline settings rather than invent them. The proof's native FFT normalization remains mandatory.

Let U_B bound lambda_max(B_0+E H_m E^T), and use the existing source-derived U_R for R. With U=U_B+U_R>0, set **tau=1e-6 U** and factor

\[
 K=B_0+E H_mE^T-\tau I.
\]

K positive definite implies lambda_min(Abar)>tau>=1e-6 lambda_max(Abar), exactly the original sufficient direction. Positive definiteness without this shift is insufficient. U_B can retain the existing inexpensive block upper bound, adjusted by max(beta,0); a negative lower bound does not invalidate that upper bound. Loose U_R can make this test inconclusive even for a well-conditioned Abar. No clipping, diagonal repair, Ewald change or threshold relaxation follows.

Numerical acceptance needs the complete factorization status, finite positive factor diagonals, unchanged symmetry check, and recorded factor/solve residual diagnostics qualified against the existing dense convention. Original dense checks include symmetry/orthogonality 1e-10 and eigen residual 1e-10*max(1,hi); a Cholesky residual has different scaling, so its arithmetic acceptance rule must be declared from backend qualification rather than relabeled as that eigen residual. Borderline pivots/large residuals return **numerically unresolved**, not passed. Do not impose an invented physical safety margin or claim an interval proof from floating pivots.

## 2. Practical backend and isolated bootstrap

Local inspection found NumPy2.4.1/SciPy1.17.0 with SuperLU, and system libsuperlu6 6.0.1. The system package is not claimed to be SciPy's embedded version. No CHOLMOD/scikit-sparse, MUMPS, PARDISO or PETSc binding appeared in the inspected active Python/system linker inventory. This is not an exhaustive search of unrelated environments.

Prefer an isolated SuiteSparse/CHOLMOD installation over a custom SPD test built around general LU. CHOLMOD supports symbolic analysis, sparse Cholesky and update/downdate; its supernodal routine computes LL^T. Critically, its C function can return true while reporting NOT_POSDEF and an incomplete factor: inspect status and minor==n, not the Boolean alone. [Official SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse), [factorization source](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/CHOLMOD/Supernodal/cholmod_super_numeric.c).

Root found `/shared/software/miniconda3/latest/condabin/mamba` and `conda`. A future bootstrap can resolve the following into a new prefix, without changing active Python:

```text
mamba create --dry-run --json --override-channels -c conda-forge \
  -p <W/J/periodic_stability_definition_v1/sparse_factorization_design_v1/environment_v1> \
  python=3.12 suitesparse scikit-sparse numpy scipy
```

This is a prospective resolution command, **not run here**. Pin its selected builds/URLs/hashes, then install only that reviewed explicit environment specification. Pin CHOLMOD, wrapper, BLAS, integer width, ordering and effective settings; force all numerical threads to one for the first prototype. The Python API changed in scikit-sparse0.5, so inspect the installed version before writing calls. Use a small C/ISO_C_BINDING interface later if needed to expose full status/symbolic memory statistics; do not assume an older Python API. [Maintainer API migration documentation](https://scikit-sparse.readthedocs.io/_/downloads/en/latest/pdf/).

SuperLU remains useful for small independent controls/linear solves, not as the chosen SPD or inertia oracle. Successful LU, an ILU preconditioner, or converged induction PCG does not establish this spectral criterion.

## 3. Work, memory and rank-three reuse

The saved pair count is 1,815,015. Keeping full 3x3 pair tensors gives at most 32,706,144 scalar entries including the diagonal. Approximate raw storage, before factor fill and temporary copies:

| Representation | Raw storage |
|---|---:|
| Full CSC, float64/int32 | 374.4 MiB |
| One symmetric triangle, float64/int32 | 187.5 MiB |
| One triangle, float64/int64 | 250.1 MiB |
| Dense float64 n by n | 9.59 GiB |

Use one triangle and streamed pair assembly; avoid Python objects for millions of entries. The native active-site map, minimum images, scaling exceptions and Ewald self term must match the inspected source. CHOLMOD symbolic analysis must report estimated factor entries/flops/front workspace and enforce a declared memory cap **before** numerical allocation. No numerical factor memory/time is known yet. Three-dimensional local graphs often have factor storage O(n^(4/3)) and work O(n²) with nested dissection, but the large neighbor count and ordering can make constants substantial; worst-case fill/storage are O(n²), work O(n³). Dense Cholesky alone has a leading term of about 1.54e13 flops here. These are scaling/planning statements, not measured costs or guarantees.

At a fixed geometry and fixed tau, a successful factor of K0=B0-tau I can be reused for H. For beta<0, three triangular solves form G=E^T K0^-1 E; K0+beta EE^T is positive definite exactly when I3+beta G is positive definite. For beta>=0, it preserves positive definiteness. Check the small matrix with the same numerical care and retain solve residuals; CHOLMOD rank updates are optional. If tau changes, diagonal shifting affects all coordinates and is not rank three. A common predeclared U across allowed H would permit reuse but may be conservative.

Changing coordinates changes real tensors and PME. Symbolic reuse is valid only while its sparsity graph remains valid; numerical factors are not transferable merely because the neighbor list is unchanged. This first design is a static feasibility check, not a claim that a global factorization every MD step is affordable. Reusing a stability decision over moving steps needs a separately justified perturbation bound or new check.

## 4. If the real part fails

B_0 may actually be indefinite while B_0+R is positive: diag(-1,2)+diag(2,0)=diag(1,2) is a simple example. A failed factor of shifted K may also reflect only the conservative shift or roundoff. Return **insufficient real-part certificate**, retain the failure details, and do not classify the full model as unstable. An instability claim would require a verified negative Rayleigh quotient of the full intended operator or an appropriately qualified full-operator calculation.

A mathematically valid later fallback is to retain a prospectively fixed set of conjugate Fourier pairs as R_S=ZZ^T, so R-R_S remains PSD. If K is nonsingular with k negative eigenvalues, a signed-pivoted LDL factor plus the small Schur matrix S=I+Z^T K^-1 Z can test K+ZZ^T. In exact arithmetic it is positive definite when S has r-k positive and k negative eigenvalues, no zeros, where r is the number of retained real mode columns. This includes cross coupling: examining only the negative-subspace diagonal is insufficient. Such a route needs a suitable pivoted symmetric-indefinite/inertia backend and its own finite numerical qualification; CHOLMOD failure alone does not supply that inertia. It is not implemented or released here.

This fallback costs r sparse solves, O(nr) storage and an r-by-r dense factor; 64 real columns alone occupy about17.5MiB. It is useful only if a modest fixed mode set resolves the bad directions. No Fourier subset, extra factorization or mode search was evaluated. A nominal saddle matrix with inverse grid convolution is not automatically sparse: the Fourier Gaussian/deconvolution inverse is nonlocal on the real grid, and zero modes require support bookkeeping. Ordinary Ritz values or PCG convergence are not interchangeable with either factorization test.

## 5. Minimal next implementation

1. Resolve/pin the isolated CHOLMOD backend and qualify its complete/error statuses on manufactured SPD, indefinite and near-threshold cases; compare against the saved small dense operators. Keep original arithmetic limits explicit, with no interval-proof prerequisite.
2. Write a standalone sparse assembler using the pinned real-pair algebra and active-coordinate map, plus an interface for the exact initialized PME upper bound. Preserve the frozen engine and all physics.
3. After a separate root release, perform one bounded symbolic/numeric feasibility attempt on an already saved box with fixed H and fixed settings; record memory, fill, work, residuals and the exact sufficient/insufficient distinction. Do not silently escalate to another method.
4. Only after success and native operator parity, integrate after current `udirect1`/direct-field PME initialization. Preserve both direct RHSs and trial channels, and invalidate on geometry/cell/grid/parameters. The old setup-time dense guard precedes PME initialization and cannot simply be replaced at the same call site. Periodic h/H forces, virial, image tails, salt composition and dynamics qualification remain separate.

**This task:** source/library inspection;13 manufactured checks, largest dimension4, one3x3 sparse factorization, 0.00844wall/0.00824CPU seconds; zero box calculations, installs, native molecular calls, fits or engine edits. The exact-rational toy verification is optional evidence only and introduces no new acceptance gate.
