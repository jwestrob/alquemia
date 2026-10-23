# One conditional joint-ceiling-removal diagnostic

**Source preparation only. No solve is released.** After the frozen component campaign is terminal and the existing full-constraint conic comparison has independently passed with Clarabel `Solved`, ask how much additional origin24 *linearized* P/B improvement is available when all four origin24 no-worsening block ceilings are removed together. Preserve the current fit, its selection, every failed origin and every acceptance rule. This counterfactual is not permission to weaken those rules, select a relaxed point, or start another fit.

Reuse the exact saved full-comparison linearized problem and immutable source/target/normalization authorities. The sole change is deleting these four Lorentz cones and their rows:

| Removed cone | Original row interval, zero-based half-open | Rows |
|---|---|---:|
| `block/old_EF` | `[0,414)` | 414 |
| `block/pair_mu` | `[414,517)` | 103 |
| `block/C_mu` | `[517,530)` | 13 |
| `block/chi` | `[530,535)` | 5 |

Retain the two P/B epigraph cones, all fifteen original force RMS cones, every original energy/force-component side, all parameter-domain/trust bounds and `0≤t≤sqrt(q_current)`. Keep the identical objective `min t`, with `q=t²`, and the original P/B denominators. The reduced problem has 17 SOCs (426 rows) plus one 856-row nonnegative cone: 1282 rows, 13 variables. It retains all 25 physical E/F gates and 847 original source inequalities, with only the four specified block inequalities removed. No per-block or per-metal sweep is allowed.

Use the exact qualified Clarabel 0.11.1 environment and all frozen parent effective settings: one thread, QDLDL, at most 200 iterations/60 seconds, original tolerances and regularization. One exclusive attempt marker precedes solver construction. Exactly one additional `solve()` call; no fallback, retry, numerical clamp or changed tolerance. The full constrained problem is never solved again.

Require the parent actual review to have status `passed_saved_conic_actual_review`, exact reviewer implementation/result/source-review pins, `certificate_checks_passed=true`, `solver_success=true` and `solver_status='Solved'`. The full result must also report `Solved`. An unsuccessful full baseline is an uninterpretable counterfactual baseline and must stop this diagnostic before an attempt. Own source review, root release and exact terminal campaign binding remain required. These pins stay unbound in source preparation.

Report raw solver status/x/s/z/gap/residuals, all retained cone and source margins, every physical gate, all six predicted block scores and changes relative to the same origin24, and the four removed-ceiling violations. Compare the full constrained and relaxed full-return q and their independently reconstructed residual norms. For the same four inherited fractions only, report exact scalar certificates under both retained and original rules; no fraction is selected or evaluated natively. Preserve the original certificate's `eligible` flag as an original-policy observation, never as a relaxed candidate authorization. No new fractions or search are introduced.

Successful floating conic solves support a local affine tradeoff diagnosis. They do not establish nonlinear improvement, a physical mechanism, global model capacity, or a rigorously rounded optimum. If the apparent gain requires larger pair/collective moment or susceptibility errors, report those increases explicitly. An unqualified relaxed return remains unqualified. Maximum native/reference/held-out calls: zero. No campaign artifacts, native source or scientific thresholds change.
