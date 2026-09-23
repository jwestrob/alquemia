# Solvent-aware response works, but does not improve this discriminator

Completed 2026-09-23, job1210579. **Close this six-source scorer branch without
expansion or refitting.** Reliable standalone solvent gradients produced valid,
affordable proposals and repaired one raw canonical ordering. They did not improve
any transferred classification over the existing three-candidate standalone pool,
and the strongest native adaptive method remains better on these same sources.
Retain the gradient/proposal interface as a research capability; no deployment or
default change is supported.

## Discriminatory utility

All six sources are available for every comparison. They are five protein groups:
the two A0AC samples are structural repeats. These consumed, deliberately difficult
PQQ examples test functional class/structural robustness, not independent direct
La/Ca affinities. No threshold was fitted to them.

| Source | Known class | Static standalone | Prior three candidates | Solvent response | New minus prior R, kcal/mol |
|---|---|---|---|---|---:|
| Q9Z4J7 | Ca | Ca | La | La | -1.573861 |
| C5AXV8 | La | La | La | La | -0.034518 |
| Q88JH5 | Ca | Ca | inconclusive | inconclusive | -1.955953 |
| A0A3F2YLY8 Ca-sample1 | La | inconclusive | La | La | +0.031006 |
| A0ACD6B9F2 Ca-sample4 | La | Ca | inconclusive | inconclusive | -1.313620 |
| A0ACD6B9F2 La-sample4 | La | Ca | La | La | -0.547697 |

The static standalone reference is valid for the static model. Its bands are
**transferred diagnostics only** for the two changed candidate models; neither has
a qualified new reference. Counts are static3correct/2wrong/1inconclusive, prior
three3/1/2, solvent response3/1/2, with zero unavailable. The native adaptive
comparator uses its own frozen canonical reference and yields5correct/1wrong;
only A0A3 Ca-sample1 is wrong. All unrounded scores and component terms are retained
in the comparison and full matrix, linked below.

The Q9(Ca)-versus-C5AX(La) ordering genuinely improves: `R(C5AX)-R(Q9)` changes
from -1.112975 to +0.426367kcal/mol. However, the full six-source class gap
`min(R_La)-max(R_Ca)` stays negative: static-9.491621, prior three-9.014910,
solvent response-8.754669kcal/mol. A0AC Ca-sample4 remains below the Ca references.
Thus the small canonical repair is insufficient evidence to calibrate or expand
this method. Larger R remains more La-like; these model contrasts are not binding
free energies or affinity probabilities.

## What ran and what the forces show

The predeclared [plan](PLAN_v2.md) and [review](PREFLIGHT_REVIEW.md) were followed.
The same source state, four physical angular modes, MACE checkpoint and standalone
xTB6.7.1 accuracy0.02 were retained. Energy and gradient were evaluated for
`MACE(vac)+GFN2(ALPB)-GFN2(vac)` with the original source/cap Jacobian. Start
SP-versus-gradient energy replay differed by at most6.43e-10kcal/mol per medium.
No native-ORCA solvent value, invented curvature or numerical DFT gradient entered.

All12 searches returned SLSQP success in3–11iterations, using3–12 distinct points
each,95points total. None hit the20iteration/40point limits. There were no SCF,
native model, analytic-gradient or geometry failures. Eleven new geometries were
retained; A0A3 Ca retained its identical starting candidate. Both metal rows scored
the common candidate pool, with exact geometry deduplication recorded.

| Source | Ca work from start | La work from start | Ca gradient norm, start→chosen | La gradient norm, start→chosen |
|---|---:|---:|---:|---:|
| Q9 | -1.888596 | -0.314735 | 23.817→0.632 | 13.886→4.834 |
| C5AX | -0.824187 | -0.789668 | 17.117→0.292 | 16.282→0.243 |
| Q88 | -2.257246 | -0.301293 | 37.119→2.109 | 33.536→29.261 |
| A0A3 Ca1 | 0 | -0.031006 | 11.668→11.668 | 4.938→0.649 |
| A0AC Ca4 | -1.571666 | -0.258046 | 31.209→0.264 | 13.664→0.612 |
| A0AC La4 | -0.914181 | -0.366483 | 27.065→0.373 | 19.334→0.273 |

Works are kcal/mol; norms are Euclidean norms of the four physical angular
**gradients**, kcal/mol/radian. Generalized force is their negative. The comparison
retains initial/chosen native, solvent and total vectors and their exact mode IDs.
Solvent-guided movement often raises native MACE energy while lowering solvent
energy more, demonstrating a distinct response rather than repetition of the
vacuum optimization. For example Q9 Ca changes by+1.345554 native and-3.234150
solvent, totaling-1.888596kcal/mol. This alone does not establish biological accuracy.

These are **finite candidates, not stationary minima**. Thirty-six intermediate
points exceeded the nonlinear displacement constraint (maximum0.803537Å), as
allowed for SLSQP trials. All chosen candidates pass the frozen final0.8Å bound
with its existing1e-7Å numerical tolerance; three touch the displacement boundary,
none the angular boundary. The final optimizer points for Q88 La and A0A3 Ca
exceeded the bound by6.27e-7 and7.91e-7Å, so the predeclared best-feasible rule
retained earlier points. This explains their large residual gradients despite
the optimizer's success message. No tolerance or stopping criterion was changed.

## Execution, cost and verification

- 228/228 new standalone calls completed:16missing baseline,190search-gradient,
 22final cross-state SP;56historical baseline cells reused.
- 94/94 new native MACE evaluations completed:83search and11cross;12start forces
 reused from36verified available source force cells. Zero new DFT calls.
- Allocation:218s×32CPU=**6976allocated core-seconds**, oneH200 reserved218s,
 200000MiB host memory requested. The execution body took212.772s.
- Actual summed call wall times:507.825s standalone (parallel),14.788s MACE.
 GPU model loading1.197s; worker span192.108s includes idle time. Peak CUDA
 allocation5313206272bytes. Local preparation/tests/reporting CPU is additional
 and unmetered. No matched-hardware speed claim follows.
- [Seven tests pass](TESTS_final_v1.txt), zero skips: actual state/mapping/cache,
 charge/recipe, gradient projection against an archived real derivative, frozen
 scope rejection, actual completed searches/matrix algebra, and an explicitly
 corrupted real-result copy checking missing-data reporting without fallback.

The full costs include all attempted calls; there were no retries. run_v2 was
preparation-only, preserved when three technical safeguards were corrected before
run_v3 launch. Original outputs, default scorer and references remain unchanged.

## Artifacts and next action

- [Unrounded comparison](../../workspaces/standalone_response_20260923/run_v3/COMPARISON.json)
- [Full common matrices and source traces](../../workspaces/standalone_response_20260923/run_v3/collection.json)
- [Actual search trajectories](../../workspaces/standalone_response_20260923/run_v3/search_collection.json)
- [Measured costs](../../workspaces/standalone_response_20260923/run_v3/COSTS.json)
- [Compact artifact pins and counts](ARTIFACTS.json), [commands](COMMANDS.md).

No further molecular work is authorized or recommended by this result. Retain the
native adaptive discriminator and reusable qualified gradient interface. A future
research question would require its own declared scope; do not widen this search,
refit on these six sources, or advertise its candidates as equilibrium structures.
