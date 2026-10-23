# Isolated CHOLMOD backend: tiny numerical qualification

The exact 56-package binary plan is installed in `/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/_envs/lanm_cholmod_v1` and passes the declared manufactured panel: **60/60 checks**, maximum dimension six. This qualifies the inspected Python API path on these artificial cases. It does not qualify any vacuum or periodic molecular matrix, native engine, large-factor memory envelope, or MD stability guard.

## Exact environment and preserved recovery

Python3.12.15 / SuiteSparse7.10.1 / libcholmod5.3.1 / scikit-sparse0.5.0 / NumPy2.5.4 / SciPy1.18.1 / pthreads OpenBLAS0.3.34. All56 installed identities, URLs, MD5/SHA256 values match the frozen explicit plan; runtime libraries and metadata are pinned in `ENVIRONMENT_AUDIT.json`. Actual OpenBLAS and OpenMP thread queries both return1. Existing environments were not changed.

The original build-glob resolution failure and too-long nested-prefix installation failure remain intact. The successful recovery changed only the prefix; it reused the same56-package plan. Numerical_v1 remains a failed author-side NumPy-Boolean JSON serialization attempt after its first tiny factorization. QUALIFY_v2 changes only scalar serialization and uses unchanged matrices/settings/criteria.

## Numerical results and failure semantics

Fixed AMD ordering and simplicial LL factorization were used. Rank modifications explicitly convert to LDL. The panel contains11 matrices, four rank3 cases, and five rejected malformed inputs:15 fresh factorizations plus3 rank mutations (19 numeric attempts including failed-v1's first factor, below the64 cap).

- Clear SPD and shifted-threshold-above cases pass. Largest checked LL reconstruction backward error is 1.4802988e-16; largest solve backward error is 7.4014683e-17, under each declared100*n*machine-epsilon bound.
- Below-threshold, singular, and indefinite LL cases raise `CholmodNotPositiveDefiniteError`, with code1/failing-column information in the message. The public `is_numeric` flag nevertheless remains true; it is not an admission criterion.
- **The indefinite rank3 downdate returns normally. Its LDL diagonal starts at -0.250006.** The independent small Schur test and explicit positive-diagonal gate reject it; the disposable factor is discarded without a solve or LL conversion. Return success must never be treated as SPD. Its explicit LDL reconstruction still closes to roundoff.
- Threshold equality and threshold+machine-epsilon remain **numerically unresolved**, including the latter's successful factorization and small backward residual. Backward accuracy does not imply forward accuracy or resolve a spectral boundary.
- The negative-Gershgorin SPD example factors successfully. The fixed B-indefinite/B+PSD-R-positive pair demonstrates why a failed real-part sufficient test cannot prove the full operator unstable.

Public scikit-sparse0.5.0 exposes no direct `minor` or `status`; there was no private-pointer access. A future native C integration must expose and check the native status/minor plus complete finite positive factors. Current Python failures preserve exception messages and public metadata.

## Meaning of the stability decision

For a supplied upper bound U, the numerical test factors K=A-1e-6*U*I. Positive K is sufficient for the inherited relative threshold when U>=lambda_max(A). These tests supply manufactured U values and check them against tiny dense references. This is ordinary double precision with explicit symmetry, eigen-residual, orthogonality, factor, and solve checks; it is not a directed-rounding certificate. No diagonal/eigenvalue repair, physical safety margin, fitted weights, or threshold changes were used.

## Recorded cost and limits

All four resolution/install attempts together used 178.690259 recorded process-wall seconds and 124.072702 child CPU seconds; failed attempts are included separately in the audit. The environment-audit body used 1.049513s wall/0.746159s CPU. Successful numerical body: 0.519241s wall/0.394159s CPU, one numerical thread; outer process wall1.012485s. Failed numerical-v1 outer wall1.057205s; its process CPU was not retained and is not inferred. These are recorded process scopes, not elapsed project time.

Zero native/PME calls, zero saved molecular matrix reads, zero scientific model/parameter changes. Sparse fill, real-space definiteness, reciprocal rescue, and actual-box memory/work remain unmeasured. A separate declared saved-matrix comparison is the next bounded numerical step; no real-box factorization is released by this result.
