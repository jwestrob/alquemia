# Existing parameters can improve the components together, but the repair is small

The saved trials do not support declaring the current water-response family
incapable of joint improvement. Of 96 previously evaluated trials, **21 reduce
both pair and collective energy errors, do not increase any of the four original
fitted blocks, and pass all 25 energy/force checks**. Twelve of those 21 also
start from an origin that passes all 25 checks. These are local improvements
relative to their own origins, not 21 independent validations or repaired models.

The remaining component errors are large. The previous fit omitted these energy
components from its objective; its failure therefore does not uniquely establish
a missing interaction term. The next useful experiment is a constrained local
step that explicitly includes them and preserves the individual force criteria.

## What was calculated

The prospectively declared [protocol](PROTOCOL.md) reuses all 673 complete
parameter vectors, 24 original Jacobian origins, 72 tangent linear programs and
96 existing trials. There were **zero new molecular evaluations, reference
calculations or parameter evaluations**. All original 531-row Jacobians replay
exactly. The five known-answer LP controls and 12 narrow arithmetic/provenance
controls pass. This is consumed-data development; no untouched validation was read.

Both bounded independent reviews pass. Singer reproduced reference joins,
component derivatives and all trial classifications. Planck independently
reconstructed the final origin, its 24 probes and four trials from 29 saved
vectors, including the three LP certificates and the numerical force-gate
metrics. Neither review repeated native calculations or solved another LP.

Two added descriptive scores are the mean squared errors of the paired P8/P9
and B8/B9 energies. P sums literal metal-water interactions; B is the remaining
collective contribution, C-P. The paired targets retain their admitted
common-parent basis and original D4 cancellation. Absolute parent total-energy
targets remain unavailable. Raw-cluster coupled-cluster energies are not
substituted for these component targets. The four old fitted blocks and all
25 old gates retain their original definitions.

## Actual finite evidence

| Existing trials | Count |
| --- | ---: |
| Both added component scores decrease | 70/96 |
| Also no increase in the original E/F block | 70/96 |
| Also no increase in all four original blocks | 55/96 |
| Also all 25 original gates pass at the trial | 21/96 |
| Also all 25 gates pass at the starting origin | 12/96 |

Only origins 3, 4, 23 and 24 pass all 25 checks. All four trials from origin 24
improve all six descriptive scores and pass the original checks. They reduce
pair-component RMS error by 0.056–0.138 kcal/mol and collective RMS error by
0.054–0.112 kcal/mol. The error remaining at the already frozen selected
candidate is about **26.56 kcal/mol for P and 22.19 kcal/mol for B**, versus
26.68 and 22.30 at the last Jacobian origin. No alternate historical candidate
has been selected from this diagnostic.

These small late improvements do not undo the earlier deterioration. At the
initial vector the respective component RMS errors were about 13.17 and 13.10
kcal/mol. A better fitted total or electrical response did not guarantee correct
energy components.

## Why the tangent direction is not a candidate

All 24 origins have a positive common-descent solution even when the four old
blocks are protected to first order. All 72 LP certificates pass; the maximum
primal-dual gap is below 7e-16. These statements concern the saved derivatives.

At the last origin, the full tangent direction has t=1.645. That is an initial
fractional slope, **not a feasible 165% error reduction**. The original E/F
weighted residual norm would rise from 0.196 to 2.179 at half step and 4.345 at
full step even in the linearized residual model. Both energy and force portions
grow. The protected directional derivatives are nearly zero, leaving positive
quadratic growth in the squared residuals. This direction must not be sent to
the molecular evaluator as an accepted repair.

The existing four final trials provide a useful smaller-step comparison: maximum
pair-component prediction errors are 0.0051–0.0206 kcal/mol, and collective
prediction errors 0.00018–0.00083 kcal/mol. This local agreement does not establish
accuracy for the much larger tangent proposal or a different origin. The last
Jacobian belongs to the penultimate accepted point, not the selected candidate.

## Decision and remaining work

Prepare a bounded proposal using the full linearized residual norms and the
individual original energy/force gates, at the last measured Jacobian origin.
Include both component errors explicitly. Verify its mathematical feasibility
before a separately declared, finite nonlinear evaluation. Keep this a local
development test; it cannot establish global representability or physical
accuracy. Preparation of that next calculation is source-only at this checkpoint.

Fresh validation, sampling-model release and protein/bulk exchange calculations
remain held. The donor-force sign disagreement is a separate unresolved problem.
All original SCF, basis, ECP and reference-method limitations remain. Quantum
collective remainders do not identify an individual classical energy term.

## Evidence and cost

[Result](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/actual_v1/RESULT.json)
contains every vector, origin, trial, LP certificate, input hash and gate status.
[Root release](ROOT_RELEASE.json) binds the exact implementation and inputs.
[Reference audit](reference_review/AUDIT.txt) preserves target compatibility and
uncertainty. The main calculation used **42.23 seconds wall, 37.64 CPU-seconds,
one thread, and zero scheduler allocation**. Source controls and subsequent saved
reviews are separate costs. No earlier chemical work was repeated.

[Reference/result review](actual_review_v1/REVIEW.json) is `cc93f38f` and used
1.827 wall/1.802 CPU-seconds. The
[independent final-origin review](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/math_review_v1/actual_v1/REVIEW.json)
is `972b73b2` and used 2.348 wall/2.322 CPU-seconds. Both bind result `ce2fd338`.
