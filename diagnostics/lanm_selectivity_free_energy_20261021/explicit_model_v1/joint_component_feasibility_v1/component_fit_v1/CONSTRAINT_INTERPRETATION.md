# What the campaign's constraints establish

**The four fixed objective ceilings are development safeguards, not physical/reference qualification thresholds.** A limited repair under those ceilings would establish only what this eight-step protocol achieved. It would not establish that the Hamiltonian cannot represent the consumed chemistry. This note reads frozen sources and the earlier saved-diagnostic report only; it reads no additional live fit vectors or outcomes and changes no gate.

## Exact ceilings and normalization

The [campaign declaration](ROOT_DECLARATION.json) fixes these **origin24** scores for every step:

| Block | Fixed ceiling | Meaning of the score |
| --- | ---: | --- |
| old_EF | 0.15340308137508737 | Original 413-residual engineering loss: equal energy, La-force, Dy-force and paired-force blocks. |
| pair_mu | 0.37816202481323913 | Mean squared Cartesian pair-dipole error / (1 D)², averaged within each of four parents and then equally across parents. |
| C_mu | 0.0052286591834414045 | Mean squared Cartesian full-shell interaction-dipole error / (1 D)² over four parents and three components. |
| chi | 0.0011670073007337461 | Mean squared error / (1 a0³)² over four state chi_zz values, using the fixed h=0.001 au energy-curvature reference. |

The original energy subblock gives equal weight to the mean squares of four balanced-exchange and four individual radial-work errors, on the 1 kcal/mol scale. Each force block averages five category-balanced rigid-force RMS errors on the 2 kcal/mol/Å scale. The two paired radial-work gates remain required even though they are not extra residual entries. Paired quantities are derived and correlated.

The parent joint objective is one quarter of the sum of the four scores. Its 531-vector concatenates half of each block's already normalized residual vector, so **each slice norm squared equals its reported block score divided by four**. CORE correctly uses fixed ceiling/4 for the four slice constraints. These are fixed-origin limits, not a monotonic per-step ratchet; individual observables may worsen while an aggregate stays below its ceiling.

None of these four numerical ceiling values is derived from an experimental error bar, a demonstrated quantum-method uncertainty or an independent validation criterion. They preserve a particular fitted compromise. Even the old_EF ceiling is additional to, and not equivalent to, passing its individual energy/force gates. The parent objective explicitly returns `response_physical_acceptance=None` and `joint_physical_acceptance=None`.

The component objective likewise uses development normalization: `q=max(J_P/711.9511599971815, J_B/497.4305961564669)`, with J_P and J_B the respective two-error mean squares at W8/W9. Those denominators are origin24 errors, not variances. The q≤0.25 milestone and both component RMS≤1 kcal/mol target are declared development goals, not new physical qualification.

## Thresholds that actually remain binding

The existing **25 model-to-provisional-reference agreement gates** are ten absolute energy errors ≤1 kcal/mol and fifteen rigid-force checks requiring both RMS≤2 and maximum body-component error≤5 kcal/mol/Å. They are genuine frozen acceptance criteria, but not confidence intervals for physical truth. They must remain separate from the fitted averages.

No physical dipole or susceptibility error tolerance was assigned. The 1 D and 1 a0³ objective scales are not such tolerances. The 0.5 a0³ state / 1 a0³ difference field checks, 1e−4 e·a0 odd-slope check, printing/CG allowances and matrix residual bounds concern numerical consistency. They cannot justify the origin24 ceilings or be substituted as relaxed physical error bars. Native finite-value, dual-response, force-consistency and normalized-SPD requirements remain unchanged.

The references retain basis, ECP/frozen-f, state and method limitations. Agreement of the consumed raw-cluster correlated contrast does not validate the pair/collective decomposition, moments or forces. There are no absolute-parent D4-complete component targets to invent. Fresh radial/donor/protein outcomes remain closed.

## One bounded next decision

First use the **already declared single saved-origin conic comparison**, unchanged, to determine whether a solver suited to this convex proposal resolves SLSQP's numerical difficulty under the same ceilings. No new solver call is authorized by this note.

If terminal synthesis still needs to distinguish ceiling restriction from local model capacity, prospectively declare **one matched saved-origin ceiling-ablation comparison**, not another fit: retain origin24, its exact 24 probes, component minimax objective/denominators, parameter/trust box, all 25 individual physical gates and numerical certification; omit only the four origin24 aggregate-score inequalities. Compare its achievable q and all four old scores with the unchanged-ceiling conic result. Do not search weights, origins, fractions or alternate tolerances. This proposed diagnostic needs its own declaration/release and does not extend the current one-solve scope.

A certified improvement available only in that counterfactual would identify a cost of the development safeguards in this local residual model. Report every lost electrical-response/old-EF score; it is not a physically accepted candidate, since electrical tolerances have not been established. Little improvement, an inconclusive solver result, or a local infeasibility certificate would still concern this Jacobian/trust region and reference set, not global Hamiltonian capacity. Any subsequent nonlinear test would be a separately bounded consumed-data decision; none follows automatically and no held-out result is needed for this diagnosis.

## Source anchors

Campaign declaration `a28ed5d3`, CORE `bcaefc3a`; parent [objective definition](../../response_matrix_export_v1/joint_repair_design_v1/objective_v1/SPEC.md) and OBJECTIVE `af7d26f6`; original EF objective `60274b54`; parent [data contract](../../response_matrix_export_v1/joint_repair_design_v1/DATA_CONTRACT.md), evaluator `51e0ba13`; earlier [saved feasibility report](../REPORT.md) and [reference audit](../reference_review/AUDIT.txt). No calculation, new fit, source/gate modification or additional outcome inspection was performed.
