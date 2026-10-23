# One component-aware bounded development fit: implementation outline

Source-only and prospective relative to the pending native half point. No fit code, objective, weights, constraints, start or budget was changed or executed. Root must declare the scientific choices below before implementation/release. This is one finite development campaign, not a repeated half-step sequence.

## Existing interfaces and the necessary distinction

Use the path aliases and exact authority table in [IMPLEMENTATION_MAP.md](IMPLEMENTATION_MAP.md).

- `D/J/evaluator_v1/EVALUATE.py` (`51e0ba13`): `Evaluator(config, cache_dir, root_release=None)`, `evaluate(values) -> {model,response,score,cells}`, `snapshot()` and `stop_event`. It supplies the qualified native physics and original scoring, with 46 metal fixtures, 10 analytic gradients, four susceptibility aliases and 48 static terms.
- `D/J/objective_v1/OBJECTIVE.py` (`af7d26f6`): delegates the original 413-entry E/F scorer and returns a 531-entry residual with four equally weighted blocks. The pair-moment, full-shell C-moment and susceptibility targets remain exact existing definitions. **P/B energies are expressly excluded.**
- `D/J/REFINE.py` (`73f36458`): hashes and preserves every vector, uses 24 central parameter probes per origin and four fixed trial slots, and selects current/trial by strict objective improvement. It reports the 25 gates but **does not require them for trial selection**. Its 24/672/24000/12 caps, source bindings, start, release and residual dimension are hard-coded.
- `D/J/DIRECTIONS.py`: scaled Gauss–Newton and BVLS enforce parameter/trust boxes only. They do not enforce physical gates or full protected-block norms. The separate failed SLSQP proposal is not a successful constrained solver that can silently replace them.

Consequently, a faithful new fit requires a versioned objective adapter and explicit selection/constraint policy, plus a versioned driver/release. Relabeling the old run or simply adding terms to its report is insufficient.

## Smallest scientific change

Keep the native evaluator's configuration, rendering, parsers, helpers and binary/cache identities unchanged as scientific authorities. Add a pure scoring layer after its existing `evaluate` output; do not put the new score into historical cached cells. From `cells[*].energy` plus the manifest's static energies, reuse the admitted literal C/P/B dictionaries and paired targets from `interaction_decomposition_v1/actual_v1/RESULT.json` (`a9ea13bf`). Every vector already supplies these energies, so this adds **zero native programs per vector**.

Preserve `old_score` exactly. Report the four new residuals in fixed order `(P8,B8,P9,B9)`, individual signed errors and separate J_P/J_B. C=P+B and W9-minus-W8 are derived reports, not additional fitted observations. No absolute-parent target, raw-cluster CC substitution, induced-only moment target or new field estimator is introduced.

A possible six-block interface is

`new_score = {objective, block_objectives, weighted_residuals, residual_layout, old_score, original_EF_checks_passed, components, constraint_metrics}`.

Root must choose the six weights and P/B engineering scales; this outline supplies no defaults. If the declared objective is a weighted sum of the existing four block scores plus J_P and J_B, with J_P=mean((e_P/s_P)^2) and similarly B, it has **535 residual entries**. For each old block with weight w, multiply its existing 531-vector slice by `2*sqrt(w)` because the old residual already carries a factor 1/2. Append `sqrt(w_P/2)*e_P/s_P` and `sqrt(w_B/2)*e_B/s_B`. Check squared-norm closure against the six block scores. Do not double-weight the old slices.

## Choices root must freeze

1. **Start:** one exact admitted vector, with its immutable origin/cell authority. A half-point success does not automatically select it for fitting. The original initial vector fails some gates; choosing it requires an explicit feasibility-restoration policy, not pretending it is a feasible constrained start. Preserve original baseline and selected-model comparisons.
2. **Objective:** scales, weights, component recovery criterion and correlation treatment. The 1 kcal/mol descriptive component scale is not an established physical tolerance. No response block has an independently justified physical error bar.
3. **Constraints and selection:** unchanged 25 gate thresholds; whether they are hard constraints on accepted iterates; and any fixed response protection budgets. Requiring every old block to be nonincreasing at every step is a modeling/optimization choice, not a physical law, and may retain the earlier fit's compromise. If constrained proposals are used, choose a solver and numerical admission prospectively; otherwise distinguish box-constrained proposals from actual gate-filtered acceptance. Report failed trials. No derivative probe may become a selected candidate.
4. **Finite algorithm:** derivative spans, trust bounds, fixed trial inventory, exact tie rule, solver failure handling, no-improvement stop and termination on substantive recovery or the hard budget. Do not inherit 24 more steps or automatic continuation by default. An iteration cap does not mean convergence.
5. **Regression:** keep the eight consumed E3/mu3 cycles as postselection diagnostic regression unless explicitly changing the training contract. Six additional trimer model calls suffice after selecting a complete 46-cell candidate; two trimer states are already present. Do not quietly train on the regression and then call it independent validation.

## Exact reuse and bounded driver changes

Make a separate driver version retaining atomic persistence, immutable per-vector artifacts, actual-source pins, no-retry failure semantics, complete denominators and stop behavior. Replace hard-coded residual dimension, release/source identifiers, initial-vector assertion, objective/constraint acceptance and finite limits only as declared. A new wrapper score must expose the new objective to the driver while keeping old score/gates independently available. The native evaluator's qualified methods remain as mapped previously.

Use a new ledger/output directory. Reuse baseline/static cells and exact already admitted per-metal cache matches without rewriting historical producers. Bind new execution policy separately from scientific cache identity. Never instantiate into the old fit cache. At an exactly reused old origin, its 24 saved derivative vectors can be rescored for the component objective without new calls; at a different start there is no corresponding saved Jacobian. Saved component derivatives cannot be transplanted to the selected final vector or the half point.

Before release, focused saved controls can establish score closure, unchanged old scores/25 gates, exact component targets/cycles, new selection and tie behavior, failure/cap accounting and initial/cache reuse. Retain all raw numerical admission and actual native qualification. A full native requalification or rerun of the previous fit is unnecessary for unchanged kernels.

## Cost and scientific stopping

One fully changed vector costs at most 56 programs; one changed-metal probe at most 28. If the original 24-probe/four-trial pattern is retained, the ceiling is **896 programs per fresh outer step**: 24×28+4×56, before cache reuse. Its former 24-step run used 21,504 programs and 11,770 seconds on 12 CPUs (141,240 allocated CPU-seconds). That is historical throughput, not a promise for a new campaign. Root must choose a finite new cap, launcher and cost review; reuse can reduce calls but cannot be assumed for unmeasured points.

Report numerical admission, original gates, absolute component errors, old response errors and consumed regression separately. A decrease from roughly 26/22 to 25/21 kcal/mol is local progress, not qualified interactions. If the declared bounded component-aware attempt cannot obtain substantive recovery without persistent gate/regression conflict, stop and examine that specific conflict before one functional change. No reference job or extra local continuation follows automatically.

Reference limits remain in [the reference audit](../reference_review/AUDIT.txt): paired common-parent electronic-plus-D4 definitions; no absolute-parent total targets; no correlated fragment/force/response qualification from the raw cluster CC contrast; retained SCF/ECP/state uncertainty. Donor-force disagreement, sampling, experimental mapping and independent validation remain separate prerequisites in the project plan.
