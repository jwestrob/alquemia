# Exact-final local constraint attribution — source/preparation only

## Question and fixed authority

Can the remaining paired P/B errors improve locally at final `feca1810ed27e323ae37e500b1d9c46f3ec2d0c34b652eeab2e2156ac62d45c8`, and which restrictions consume that flexibility: parameter-box faces, the original 25 energy/force gates, or preservation of the current electrical block norms?

Root's [scientific decision](../terminal_synthesis_v2/ROOT_DECISION.json), SHA256 `5ac3f47242acd8738f4a1d1e34fd90043b249675ece785e68dd2f8800a9b4941`, authorizes this source preparation and metadata inventory. **It authorizes no native probes or numerical solves.** All execution/release fields remain null. The original eight-step campaign and all failures stay immutable.

[INPUTS.json](INPUTS.json) binds admitted terminal result `4a374ff3`, review `13f23ffd`, selected artifact `d7ce78d0`, the exact old evaluator/manifest/CORE/domain/targets, reuse registry, last proposal and Clarabel settings/environment. The synthesis's earlier suggestion of signed electrical nonincrease rows is refined here to **three aggregate norm constraints**, avoiding per-coordinate restrictions on symmetry noise. Signed electrical predictions remain required outputs.

## Finite derivative inventory

Use the inherited `CORE.probe_points` algorithm, literally: for each of twelve parameters in declared order, minus then plus 0.005 in normalized box coordinates, clipped to [0,1]. Preserve its binary64 arithmetic and hash the complete twelve-value vector. Divide response differences by the **actual** normalized endpoint span. A clipped pair is an asymmetric secant, not a symmetric derivative. In particular the negative Dy_v probe reaches −1; its positive endpoint is separated by the full nominal increment. No extra stencil, extrapolation or changed box is implied.

`PREPARE.py` executes only AST-extracted `require/token/checked/probe_points` from frozen CORE and `require/checked/effective` from the evaluator. It does not import either execution module, instantiate a solver, invoke native code, read old mutable ledgers, or replay derived LAPACK quantities. It joins the admitted immutable historical registry with terminal metadata and checks complete-vector keys plus exact per-metal effective `[a,b,c,d,R,epsilon]` keys. The selected body is the preferred unchanged-metal authority. A future executor must still verify the matched immutable cell/cache/render provenance before import.

Maximum scope is 24 derivative vectors × 56 programs = **1,344 programs before reuse**. Each vector has the original 46 cells: 23 per metal, with five gradient programs plus 23 analyze programs per metal. Exact unchanged-metal and any complete-vector reuse reduce the prospective count; the inventory records that reduction. No new geometries, fields, quantum calls, trimer targets, heldouts or fresh radial/donor/protein outcomes are included.

The completed metadata preparation finds **0/24 complete-vector matches**, twelve new effective La states and twelve new effective Dy states. Every unchanged-metal half matches the admitted selected state. The prospective total is therefore **672 new programs** (24×28), subject to exact cell import verification before any future execution. Dy_v's minus/plus normalized distances are 0.000170593745372649/0.005000000000000004, with total span 0.005170593745372654; every other probe pair is unclipped. Preparation consumed 0.603524 s wall/0.584238 CPU-s on one thread and zero native/solver/reference calls. Its immutable inventory hash is `baeaefc586cc876990cc99e39c02795b831f166ba4a493a6a28f1446a91a1542`.

The last saved Jacobian in `PROPOSAL_08` is at `a879be2f…`, the selected step7 point. Its saved trials may diagnose off-center prediction quality and restriction patterns; it is **not** the final-point Jacobian. Missing final derivatives remain missing until a separate release and actual admission.

## One linearization, four fixed comparisons

After all required final derivatives are actually admitted, construct one common linearization using frozen CORE physical413/weighted531 and literal paired P8/B8/P9/B9 definitions. Do not change the common-basis electronic-plus-recorded-D4 references or invent absolute-parent targets. The original two component squared-error denominators stay `J_P0=711.9511599971815` and `J_B0=497.4305961564669`.

Let `z` be twelve normalized parameter increments, `t≥0`, and `q=t²`. For the two component vectors, minimize t subject to

`||e_P + J_P z||₂ ≤ sqrt(2 J_P0) t`,

`||e_B + J_B z||₂ ≤ sqrt(2 J_B0) t`.

Every case retains the **same local trust region |z_j|≤0.025**, the same derivative/scoring data, two component cones, and objective. The raw epigraph and original q are retained. For numerical reporting only, construct `q_aux=max(recomputed P/B squared-error ratios)+1e−8` at the unchanged physical vector, using the previously declared auxiliary-only convention. The common upper objective bound is the center's similarly constructed q; it introduces no physical tolerance or parameter change. Report any raw epigraph violation before the auxiliary reconstruction.

| Case, fixed execution order | Parameter box | Original 25 physical gates | Current electrical norms |
|---|---|---|---|
| **A — constrained mechanics** | retain | retain | report |
| **B — electrical tradeoff** | retain | retain | retain pair, C and χ separately |
| **C — box-face sensitivity** | remove only box faces; trust remains | retain | report |
| **D — physical-gate sensitivity** | retain | remove only the 25 gates | report |

A→B measures the restriction introduced by preserving present electrical performance. C versus A attributes sensitivity to box faces; D versus A attributes sensitivity to original physical gates. These are controlled local comparisons, not a unique additive decomposition of causal importance. **C/D predictions are diagnostic only**; outside-box or failed-physical vectors cannot become candidate evaluations or recommendations for execution.

### Exact retained physical and electrical rows

The original ten energy errors retain absolute limits 1 kcal/mol, including both paired radial controls. Each of fifteen original body-force groups retains its frozen weighted RMS ≤2 kcal/mol/Å and every projected component within ±5. Use the same rigid-body/torque weights and source ordering; do not substitute Cartesian-force RMS or only the aggregate E/F block. Thus A/B/C retain the inherited fifteen force SOCs and their component inequalities, plus energy inequalities. D removes those physical rows as a group and reports their predicted violations.

For B, use the original weighted531 slices for pair/C/χ, with their original factors. Require each predicted slice norm at most its **final-point** slice norm. Equivalently protect the three existing block scores individually, with no old-E/F norm ceiling and no origin24 electrical ceiling. This adds exactly three SOCs. These are Pareto-development comparisons, **not experimental electrical accuracy tolerances**. Do not protect each near-zero Cartesian component separately. The physical gate counts/limits stay unchanged.

Box faces are those of the frozen `PARAMETER_DOMAIN.json`. La_u and Dy_v are the near-edge box coordinates; map native coefficients using c=S(u+v)/4 and d=S(u−v)/4. C's algebraic removal of those faces is not a native domain expansion or a statement about physical saturation.

## Numerical settings, output and stop rules

Prospective solve cap: **four independent one-attempt solves**, A/B/C/D in that order, one numerical thread, with the existing pinned Clarabel 0.11.1 environment and every setting from `SETTINGS.json` (`max_iter=200`, `time_limit=60 s`, `max_threads=1`, original tolerances/regularization). No warm-start trajectory, fit selection, retry, fallback, tolerance change, solver installation or follow-up native evaluation is part of this panel. This preparation performs zero solves.

Keep full raw solver status, primal/dual variables, residuals, gaps, iteration/cost receipts and matrices/row-family identities. Independently recompute affine residual norms and scalar constraint margins using the established compensated-sum certificate convention; report exact-sign flags separately from numerical feasibility at 1e−8. Gap ≤1e−7 is a reported numerical diagnostic, not a replacement for raw status or an optimality theorem. Retain failed attempts and do not infer infeasibility from `AlmostSolved`, solver failure or no measured gain. Inspect dual face/row sensitivities only with their stationarity/complementarity and residual limitations. Never relabel an auxiliary epigraph repair as solver convergence.

Required per-case outputs: P8/B8/P9/B9 and P/B RMS/q; all 25 predicted physical metrics; all four original block scores; full signed four χ states, their two trimer-minus-dimer increments and Dy−La increment difference; all four C dipole vectors and W9 Dy−La contrast; pair/C/χ RMS in original units; all twelve proposed coordinates, box/trust distances, mapped native coefficients; original/raw and auxiliary q certificates; and full/half-step predictions at the same direction (arithmetic only, no native selection). Preserve the center as a feasible reference point in each retained constraint family and report its margins.

The panel distinguishes restrictions in one measured local affine model. A lack of feasible improvement may reflect the local span, clamped secant quality or numerical/trust-region limitations and is not a global Hamiltonian no-go. A feasible diagnostic direction remains a prediction requiring a separately declared nonlinear check. Nothing here authorizes another eight-step campaign, new interaction terms, broader bounds or production/free-energy sampling.
