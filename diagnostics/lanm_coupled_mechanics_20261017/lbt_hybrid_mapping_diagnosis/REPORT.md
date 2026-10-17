# Hybrid origin-map failure: derived-coefficient roundoff

Actual failed-host diagnostic12984 completes all four origins, two repeats each, with zero molecular calls/model loads. All coordinates, identities, charges and types match exactly. Only the baseline7CCN first cap differs:16 Jacobian coefficients duplicated in the cap/core views,32 scalar leaves total, maximum1.1102230246251565e-16 and4 binary64 ULPs. All repeated calculations are exact; the other three maps match completely.

The same cap norm differs by one ULP between NumPy dot/norm and sum-of-squares arithmetic on this host. This supports a floating-point recomputation explanation; it does not identify the exact hardware instruction responsible or establish a universal cross-host error bound. The original passing CPUpreflight12972 and failed H200scout12975 remain immutable.

A versioned recovery will permit only a small declared roundoff envelope in derived cap Jacobians. Coordinates/taskXYZ and all other leaves remain exact; physical, state, energy and finite-difference gates stay unchanged. Actual-host zero-call replay and meaningful negative controls must pass before the unchanged nine-cell molecular scout. No hybrid energy or derivative result is established by this diagnostic.

Actual cost:1 scheduler second,12 allocated CPUs, no GPU; application0.548372s. The12CPU request matches the original launch environment. GFN backend import then sets numerical thread variables to1, exactly as in the failed evaluator. Runtime/import provenance and full differences are preserved in result/.
