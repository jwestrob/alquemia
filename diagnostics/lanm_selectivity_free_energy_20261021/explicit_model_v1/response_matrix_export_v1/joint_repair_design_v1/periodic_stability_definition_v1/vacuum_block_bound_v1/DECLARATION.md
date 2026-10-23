# Fixed vacuum block bounds

## Question and scope

Do simple block-Gershgorin bounds establish the existing normalized positive-definiteness criterion on the six unique nonzero-response zero-field matrices in the admitted fresh 57-program CPU qualification? Use only `response_matrix_001.dat` for the La/Dy dimer, near two-water fragment, and W9 center. Sequence 2 is excluded. No native calls, fit, new coefficients, partition search, weight search, periodic calculation or engine change is authorized by this diagnostic.

This declaration and exact input/source pins precede numerical evaluation. Run once with one numerical thread. Retain every failure. Actual wall and CPU time are recorded. The six combined profiles, input coordinates, native polarizabilities and matrix mappings remain frozen.

## Two fixed partitions

1. Each native polarizable site is one block of its three Cartesian coordinates.
2. Each intact water molecule is one block of its nine coordinates; the single metal is a separate three-coordinate block. Recover water membership from the admitted literal XYZ bond graph and verify types 36/37, two hydrogens per oxygen, reciprocal bonds and complete disjoint coverage. Do not cluster waters or optimize block membership.

## Calculation and unchanged criterion

Let `B = sqrt(alpha) A sqrt(alpha)`, repeating each admitted site's alpha for its three coordinates. Record the maximum antisymmetry and its spectral norm. For the floating-point block calculation use `S = (B + B.T)/2` and report its dense eigenvalues beside the saved normalized dense spectrum. This explicitly exposes the minute numerical symmetrization; it is not a new Hamiltonian or a rigorous interval enclosure of the raw matrix.

For each partition block i, calculate its diagonal-block eigenvalue extremes and `R_i = sum(j != i) ||S_ij||_2`, where each off-diagonal spectral norm is computed by SVD. Record all radii and off-block norms, then

`L = min_i (lambda_min(S_ii) - R_i)`

`U = max_i (lambda_max(S_ii) + R_i)`.

For an exactly symmetric matrix, these bounds follow by applying the block eigenvector equation to the block with largest Euclidean norm. They enclose every eigenvalue, without a sign assumption on individual couplings. Accept the *descriptive bound criterion* only when **L > 1e-6 U**, with positive U. The threshold is identical to the frozen native normalized-spectrum criterion. No clipping, epsilon allowance or adjustment enters this criterion. Compare L/U to the actual saved minimum/maximum eigenvalue ratio and show lower/upper gaps. The inherited scout arithmetic rule may check saved-versus-recomputed dense roundoff and floating enclosure; it never changes L, U or the criterion.

These are exploratory floating-point evaluations, not directed-rounding certificates, proofs of periodic stability, physical validation, or a replacement for the admitted dense native test. A negative lower bound means this sufficient bound fails, not that the matrix is unstable. Report atom/molecule usefulness separately for one, two and nine waters and both metals.
