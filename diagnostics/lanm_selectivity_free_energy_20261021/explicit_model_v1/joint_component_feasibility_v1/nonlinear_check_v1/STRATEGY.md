# What the pending half-point pilot can decide

Prospective scientific strategy, 22 October 2026. The native half-point result has not been inspected for this note. No new fit, parameter choice, reference calculation or continuation is released here.

## Evidence and purpose

The [saved feasibility diagnostic](../REPORT.md) contains 21 existing trials that improve both paired component scores, do not increase the four old objective blocks and pass all 25 E/F gates; 12 also start at gate-passing origins. Thus there is some local flexibility. Nevertheless, fitting the old objective worsened pair/collective RMS errors from about **13.17/13.10 to 26.56/22.19 kcal/mol**. These component energies were absent from that objective.

The [constrained proposal report](../constrained_step_v1/REPORT.md) preserves the SLSQP failure. Its fixed half point is separately feasible in the saved linearized problem: P/B RMS errors are predicted to decrease from **26.68/22.30 to 25.59/21.39 kcal/mol**, about 4%. The source result is `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/constrained_step_v1/actual_v1/RESULT.json`, SHA-256 `0338640df3045303b959f22a3e729cffce7fd5786db8e1970cb0d442096b5d5d`.

The native pilot tests whether that one prediction translates into finite component improvement while preserving the original gates and response blocks. It measures local prediction error. It cannot establish accurate interactions, cure the trimer regression, establish a successful optimizer or qualify a sampling model. Even perfect agreement leaves large errors. Some predicted margins are small: the pair-moment residual norm improves by only 1.99e-7.

## What remains unidentified

Objective omission is established; fundamental lack of model flexibility is not. Local improvements neither prove that reweighting can reach acceptable accuracy nor support global impossibility. Quantum B includes all nonpair orders and is not a classical polarization term. No unique charge-transfer, Pauli or stiffness explanation follows.

The paired common-parent electronic-plus-D4 definitions must remain intact. Absolute parent total-energy targets are unavailable, and raw cluster CC energies cannot substitute for CP/P/B or force/response targets. Reference SCF, basis/ECP/frozen-f and electronic-state limits remain. The unresolved donor mechanics are separate from this water pilot.

## Finite decision tree

1. **Native admission and the declared finite-improvement conditions pass:** close the local feasibility sequence. Prepare one prospectively bounded component-aware development fit, with explicit P/B targets, E/F and response terms, fixed scales, cost cap and stopping rule. Compare the result with both preserved starting and selected models. A few more percent improvement alone is not interaction qualification.
2. **Native admission passes but the prediction or a protected condition fails:** identify the failed observable and prediction discrepancy. Do not automatically halve again or change a tolerance. Use this evidence when deciding the one bounded development fit; failure of this point does not establish global family failure.
3. **Native admission fails:** retain this point as numerically unavailable. Diagnose the specific failure before further model work. Do not infer physical inadequacy from an implementation/solve failure alone.
4. **A properly component-aware bounded fit still leaves large errors or worsens consumed trimer behavior:** stop small continuations. Identify the persistent force/component/response conflict and use it to justify one concrete functional change. Do not add terms merely because the earlier objective omitted targets.

Before a fit, root must define what constitutes substantive component recovery and an eventual acceptable interaction error, using the intended exchange resolution and reference limitations. Engineering scales are not physical uncertainty bounds. Keep the original 25 gates unchanged; separately declare any response protection budgets. Preserve consumed trimer E3/mu3 as the diagnostic regression unless a new training contract explicitly changes that role. Fresh outcomes remain closed.

The highest-value next development decision is whether a complete component-aware objective can obtain materially better interactions within the present family. This follows the interaction-qualification stage of the project [PLAN](../../../PLAN.md); it does not justify bulk/protein production or repeated slow reference campaigns.
