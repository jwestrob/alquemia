# Fixed comparison-matrix follow-up

The first unweighted-bound attempt remains unchanged. On its same six sequence-001 vacuum snapshots and exactly its two fixed partitions, test whether positive block weights can remove the observed conservatism. Reuse the pinned diagonal-block eigenvalue extremes and off-block spectral norms from that saved calculation; do not recompute or optimize partitions, block norms, native inputs or coefficients. Verify the original raw snapshot and geometry pins. One thread, zero native calls, no optimizer.

For each partition construct the symmetric comparison matrix C with `C_ii = lambda_min(S_ii)` and `C_ij = -||S_ij||_2`. Here S is the same normalized symmetric matrix defined in the original declaration. Calculate C's small dense spectrum once. For an exact symmetric matrix, writing y_i = ||x_i|| gives `x.T S x >= y.T C y`, hence `lambda_min(S) >= lambda_min(C)`.

If `lambda_min(C) > 0`, solve `C w = 1` once and require all finite `w_i > 0`. A positive-definite symmetric matrix with nonpositive off-diagonal entries has a nonnegative inverse; this construction supplies positive weights without optimization. Report the solve residual, every weight, and

`L_weighted = min_i ((C w)_i / w_i)`

`U_weighted = max_i (lambda_max(S_ii) + sum(j != i) ||S_ij||_2 w_j / w_i)`.

Apply exactly the frozen criterion `L_weighted > 1e-6 U_weighted` with positive U. No threshold adjustment, clipping, retry, alternative right side or weight refinement is allowed.

If `lambda_min(C) <= 0`, do not solve for weights. Report that no positive-weight block-Gershgorin witness can establish a positive lower bound using this exact comparison matrix: its minimum eigenvalue has a nonnegative eigenvector v, so positive Cw and w would contradict `v.T Cw = lambda_min(C) v.T w <= 0`. Thus failure here separates loss from scalar block-norm compression from a merely poor choice of unweighted row sums. It does not imply instability of S.

All numerical results remain ordinary floating-point exploration. Neither the spectral sign nor the weight solve is a directed-rounding certificate or an admitted native guard. Compare against the original saved dense spectrum and retain every failure. The only arithmetic allowance is the already inherited floating-point check of numerical enclosure; it does not alter C, its sign, weights or the 1e-6 criterion. All twelve comparison spectra and at most twelve conditional linear solves are finite and predeclared.
