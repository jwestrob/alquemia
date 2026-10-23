# Local feasibility of all twenty-five unchanged gates

**No point in the declared local linear model and ±0.025 trust box can satisfy all twenty-five gates.** Even the necessary energy/component-only relaxation has optimal threshold multiplier **3.213222329**, where the original thresholds require at most1. This conclusion concerns the saved selected-point linearization, not the nonlinear physical family.

## Numerical certificate

HiGHS dual simplex terminates optimal after 3 iterations. Primal violation is 0; dual stationarity residual 7.53e-15; complementarity residual 5.41e-17; primal–dual gap -4.44e-16. The sparse dual certificate independently reproduces the lower bound3.2132223286572854 without another solver call.

The active scientific constraints are the central-dimer and W9 aquo8-referenced exchange errors. At the minimax point both remain−3.213222329kcal/mol. Five of the six displacement variables are on trust-box faces; original parameter bounds are included, and no claim extends beyond this frozen neighborhood.

All44 source/input/weight reconstruction checks pass: eight raw energy rows, two paired-work differences,405 raw body components and fifteen category-balanced RMS values. Exact scientific thresholds remain energy1, RMS2, component5. No threshold was changed.

## Predicted minimax point — not evaluated by a molecular program

Normalized parameter displacement: `[0.012605058049302691, -0.025, 0.025, 0.025, -0.025, 0.025]`.

Predicted physical parameter vector in the frozen six-variable order: `[3.9333809051695616, 0.1653770415255594, 3.6650821657333443, 0.16582541153280464, 0.15253146405473372, 0.2118129058341702]`.

This point passes **6/10 energy** and **9/15 force** gates. Its projected weighted loss is1.888143976; minimizing the largest gate violation is distinct from minimizing that loss.

| Energy gate | Predicted error (kcal/mol) | Pass |
|---|---:|---|
| balanced_exchange/water_minus | -2.733627 | False |
| balanced_exchange/water_central | -3.213222 | False |
| balanced_exchange/water_plus | -2.664275 | False |
| balanced_exchange/aquo9_central | -3.213222 | False |
| individual_radial_work/La_minus | 0.091526 | True |
| individual_radial_work/La_plus | -0.642825 | True |
| individual_radial_work/Dy_minus | 0.571121 | True |
| individual_radial_work/Dy_plus | -0.093877 | True |
| paired_radial_work/minus | 0.479595 | True |
| paired_radial_work/plus | 0.548948 | True |

| Force gate | Predicted RMS | Predicted maximum component | Pass |
|---|---:|---:|---|
| La/water_minus | 0.559633 | 1.069532 | True |
| La/water_central | 1.463974 | 2.906575 | True |
| La/water_plus | 2.581314 | 5.115842 | False |
| La/aquo8_central | 1.127243 | 3.430926 | True |
| La/aquo9_central | 2.160447 | 4.082768 | False |
| Dy/water_minus | 3.463886 | 6.872337 | False |
| Dy/water_central | 0.801225 | 1.590489 | True |
| Dy/water_plus | 0.030616 | 0.059956 | True |
| Dy/aquo8_central | 0.316840 | 0.708838 | True |
| Dy/aquo9_central | 2.759078 | 6.005739 | False |
| paired/water_minus | 4.012281 | 7.941869 | False |
| paired/water_central | 0.662754 | 1.316086 | True |
| paired/water_plus | 2.550917 | 5.055886 | False |
| paired/aquo8_central | 1.096213 | 3.425571 | True |
| paired/aquo9_central | 1.204056 | 2.849202 | True |

SLSQP was not run: adding RMS constraints cannot restore feasibility when the necessary LP is already certified infeasible. Every RMS ratio at this LP point is also below its3.21322 common multiplier, so the point attains that same relaxed multiplier for the full local gate set. It does not satisfy the original gates.

## Limits and provenance

The Jacobian uses normalized parameter halfsteps0.005. Linear predictions at a0.025 trust displacement are not evaluated nonlinear responses. This local certificate neither excludes a distant nonlinear solution nor establishes which physical interaction is missing. W8/W9 geometry and coordination changes remain confounded. The controls are correlated engineering criteria, not a statistical likelihood.

All completed fit/raw-review parent pins remain unchanged. The executed source and RESULT are immutable. An independent source reviewer confirmed weighting and dual conventions, and suggested explicitly requiring in-box status for a future feasible classification. That guard is in unexecuted SOLVE_v2.py; the executed certified-infeasible branch is unaffected and was not rerun.

The sparse certificate exporter had one preserved syntax failure before calculations, then completed; it made no solver call. Main diagnostic and certificate arithmetic used one thread with HiGHS threads=1 and BLAS/OMP1: 7.343s wall / 7.324CPU s, plus0.045s reported shell time for the pre-calculation syntax failure. Zero native calls, new reference outcomes, candidate changes or automatic evaluations.

Exact artifacts: RESULT.json SHA256 `3de565b490a5a8337f6860775bbb64b26dc18eb469b161e600345121d6bb3554`; SPARSE_DUAL_CERTIFICATE.json SHA256 `873d7e57708e6586d575eed54d95b5c7fe50c76cdbf3ba46728b8cfcef63795c`.
