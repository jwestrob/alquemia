# Saved geometric-tail result review

**Passed for the two released geometry-only evaluations.** Both calculator processes exited zero with empty stderr, used the reviewed implementation and exact declared configurations, and retained all scope limitations. The saved La and Dy geometric bound records are identical. This review checks saved receipts and output structure; it does not recompute a tail or rerun either process.

For the fixed 49.324 Å cubic boxes, 3,980 primary waters and inherited literal metal–oxygen distance bound D = 43.492549421087 Å, the reported omitted-tail upper bounds are:

| Retained image cube | T0, dimensionless (approximately) | T1, Å (approximately) |
| --- | ---: | ---: |
| N = 1 | 1.11155965 × 10⁻⁸⁰ | 2.38269366 × 10⁻⁷⁸ |
| N = 2 | 2.56200665 × 10⁻¹⁶⁵ | 7.68057753 × 10⁻¹⁶³ |
| N = 3 | 2.94922009 × 10⁻²⁵⁰ | 1.13609634 × 10⁻²⁴⁷ |

Exact outward endpoints remain in the pinned outputs. N = 1 retains every image with integer index max-norm at most one, including adjacent images; these results do not justify dropping those images or using a minimum-image-only definition. The very small distant-image bounds concern the specified radial environment sum and its distance-weighted counterpart at these fixed geometries. No acceptance threshold or new cutoff is selected.

D comes from the pinned prior representative certificate. The root release records successful nested source/configuration/geometry/inventory hash verification. Neither the calculator nor this result review independently recomputes D from coordinates, and no future-configuration domain is admitted.

No response coefficients, retained sums, periodic dipoles or inverse-response bound were supplied. Therefore h/H errors, total self-consistent force errors, stability and a periodic Hamiltonian remain unqualified. Tiny geometric tails alone establish none of those properties.

The original terminal receipt records two serial calculator programs, 0.614873487 s wall time and 0.603294 s child CPU time, with zero native calls and zero scheduler-allocated CPU seconds. The review adds only saved-file checks. `REVIEW.json` pins both outcomes, the terminal receipt, root release, source review and configurations.
