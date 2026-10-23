# One static supernodal factor attempt — source design

## Decision and scope

Prepare one **fixed supernodal LL factorization with AMD ordering** in the existing pinned CHOLMOD environment. This is a practical next numerical implementation to qualify, not a released calculation. No factorization, new graph analysis, solver search, installation, native call, candidate composition or engine change was performed for this design.

The existing one-box symbolic result remains a resource rejection under its original 4 GiB planning budget. It is not an instability result. Any later attempt requires the fit to be terminal, the exact reviewed candidate and BIO composition to be bound, and a bound from the initialized native PME state. It would test one static operator; it would not establish factorization at every MD step as practical.

## Why change the numerical mode

The installed `cholmod.h` describes simplicial `rowfac` as incremental factorization with elimination-tree dependencies. Its supernodal module explicitly uses BLAS/LAPACK, including `dgemm`, `dsyrk`, `dtrsm` and `dpotrf`; `cholmod_super_numeric` is documented as a BLAS-based LL method. The installed library imports those symbols and links the pinned pthreads OpenBLAS and libgomp libraries. This establishes an available blocked computation path, **not measured parallel speedup**. The available local header and wrapper declarations do not establish useful BLAS threading for simplicial factorization, or prove every simplicial operation is serial. Full numerical-kernel source is not installed.

Use `CholeskyFactor(K_lower, lower=True, order='amd', sym_kind='sym', supernodal_mode='supernodal')`, followed by exactly one `factorize(K_lower, ldl=False, beta=0.)`. The diagonal shift is already present in `K_lower`. Preserve the installed defaults for other settings. Record their pinned library/wrapper provenance and available runtime state; do not silently choose `auto`, another ordering, LDL, diagonal bounds, repair, or an alternative factor after failure. Supernodal symbolic preparation is part of the eventual bounded attempt; the prior simplicial analysis is not a supernodal allocation forecast.

## Storage, work and proposed resource envelope

The admitted graph has 35,874 coordinates, 16,371,009 lower-CSC entries and **399,831,681 simplicial symbolic L entries**. Input CSC arrays occupy 196,595,608 bytes. The packed simplicial factor-array estimate is 4,798,841,168 bytes (4.469 GiB); the earlier forecast with explicit extra-storage allowances is 10,264,187,872 bytes (9.559 GiB). Neither is measured supernodal peak memory. Relaxed supernode amalgamation, frontal workspace, factor conversion and a copied sparse L can increase storage. The public wrapper does not expose `xsize`, `ssize`, `maxcsize` or Common peak-memory counters.

The saved sum of squared column counts is 5,942,759,498,717. This is a structural work proxy, not observed CHOLMOD `fl`, elapsed time or an accuracy estimate. Sparse factorization is fill-dependent; at this graph's approximately 62% of dense lower-triangle occupancy it must not be described as near-linear. The diagnostic phase uses O(nnz(K)+nnz(L)+n) storage and O(k[ nnz(K)+nnz(L) ]) sparse work for fixed k=10 probes, without a dense n-by-n matrix.

Propose one worker with **32 GiB hard/soft address-space limit, 900 s wall limit, and 10,800 aggregate process CPU-seconds**, including imports, assembly, supernodal symbolic/numeric work, extraction and checks. Use 12 pthreads OpenBLAS threads; set OpenMP and other numerical pools to one. Verify actual OpenBLAS thread count 12 and OpenMP maximum 1 before work; record host, PID, limits, actual wall/CPU, peak RSS and phase. A separate nonnumerical supervisor kills only its own process group on timeout. Start only after root confirms the fit is terminal and the shared-host envelope is available. The previous memory snapshot showed ample host memory, but it is not a reservation; recheck at release/runtime.

The 32 GiB proposal allows substantially more than the existing 9.559 GiB forecast, including one copied sparse L and potential conversion temporaries. It is an **engineering budget, not a proven upper bound or completion promise**. Memory/timeout/extraction failure remains a preserved resource failure or inconclusive check. There is no automatic budget increase or retry. Changing to supernodal does not itself release this envelope.

## Unchanged stability question and required bindings

On the actual positive-alpha induced degrees of freedom, write the normalized operator as `Abar = B0 + H + R`, where the proved initialized PME contribution `R` is positive semidefinite. Obtain valid upper bounds `U_B` for the actual candidate `B0+H`, and `U_R` from the exact initialized PME weights/grid. Require finite `U=U_B+U_R>0`, and retain **tau = 1e-6 U**. Factor `K=B0+H-tau I` once. Positive K is a sufficient route to the inherited strict relative threshold for Abar. No eigenvalue/diagonal repair or extra physical safety margin is introduced.

Failure to establish positive K does **not** establish that Abar fails: the reciprocal term can stabilize an indefinite real-space part. Failure must therefore be reported as failure of this sufficient split test, unless independent full-operator evidence says more. Successful factors and backward checks are ordinary-double numerical evidence, not interval certificates or forward-error guarantees near the threshold. A warning or unresolved numerical boundary is not promoted to a pass; raw `rcond` is a factor-diagonal estimate, not the eigenvalue ratio.

Required before a runnable configuration: terminal candidate result/config/hash and review; explicit BIO composition and exact alpha/damping/H; geometry, cell, cutoff/split and active-site map; initialized PME export with its source review, nonnegative multiplier/adjoint conventions and finite upper bound; fixed mode qualification and root release. The saved graph used literal old BIO parameters, H=0 and diagnostic split/cutoff. Those values cannot substitute for candidate/PME binding. The local metal H changes the diagonal and h changes the RHS; neither licenses omission of its exact candidate values or the native parity checks owned separately.

## Residual checks without forming L L-transpose

After an exception-free LL factorization, require the expected supernodal/LL flags and finite positive extracted diagonal; `is_numeric` alone is insufficient. Obtain one documented copied CSC factor using `get_factor(kind='LL', lower=True)` and permutation `p=get_perm()`, where `P=I[p]`. Its copy/conversion storage stays inside the hard cap. Verify the original factor remains suitable for solving in the tiny controls. Never materialize `L @ L.T`, K dense, or a dense permutation matrix.

For a vector v, apply K from its lower CSC as `K_lower@v + K_lower.T@v - diag(K)*v`. Compute its infinity norm from sparse absolute row/column sums, subtracting the diagonal once. With `vp=v[p]`, compare `L@(L.T@vp)` to `(K@v)[p]`. Require

`||factor_action - permuted_K_action||inf / (max(1,||K||inf) ||v||inf) <= 100*n*eps`.

For the same ten vectors as right-hand sides, call the original factor's `solve(b, system='A')` and require

`||K*x-b||inf / max(tiny, ||K||inf ||x||inf + ||b||inf) <= 100*n*eps`.

Freeze ten probes before the attempt: three uniform Cartesian-coordinate vectors in normalized DOFs, three basis vectors on the metal coordinates, and four Rademacher vectors from explicitly constructed NumPy `Generator(PCG64(20261022))` with its pinned version. Process one vector at a time. All intermediates and reported norms must be finite; record every residual and maximum, not only a pass flag. Preserve the inherited symmetry checks at operator construction/admission: CHOLMOD consumes one triangle and does not check symmetry itself.

These are **sampled factor-action residuals and solve backward errors**. They do not equal the tiny-panel whole-factor matrix-norm residual and cannot certify an unsampled operator-error norm. State this limitation explicitly. They supplement a standard trusted factorization path; no mandatory interval/inverse-verification framework is proposed. The known indefinite-downdate return-success example remains preserved; this attempt performs no updates/downdates and never accepts success flags alone.

## Minimum finite controls before the mode change

Reuse the exact eleven manufactured matrices in the existing numerical declaration, maximum dimension six, with fixed supernodal LL and the unchanged symmetry, shift, error and boundary criteria. No real box or new directions are needed. Confirm expected SPD/nonpositive behavior and preserve exact/roundoff-boundary cases as unresolved. Keep the existing malformed-input rejection checks; no replacement scientific threshold.

On these tiny factors, compare the new sparse action residual calculation to a direct dense LL-transpose check; verify permutation, copied-L independence, finite positive diagonals, original factor state and solve readback. One deliberate 10% change to a copied positive diagonal in the diagonal-SPD fixture must be rejected by the action check; do not alter the native factor. No new rank-update tests are needed because the proposed static path does not use updates. This finite mode/API check may run with one thread after declaration/release; it does not benchmark 12-thread performance. The future actual worker must independently read back its declared thread settings and perform all unchanged actual residual checks.

**Status:** source design only. No supernodal controls, new symbolic analysis, numeric factorization, solve, native call or actual candidate binding has occurred in this task. Exact authorities and installed API evidence are adjacent in `SOURCES.json` and `API_EVIDENCE.txt`.
