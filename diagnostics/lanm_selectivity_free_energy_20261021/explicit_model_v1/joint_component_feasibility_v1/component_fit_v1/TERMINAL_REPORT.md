# Component-aware fit: modest recovery, interactions remain inaccurate

The bounded fit completed all eight declared steps. Pair/collective RMS errors fell from **26.68/22.30 to 22.81/19.01 kcal/mol**, while all 25 original physical checks and four fixed development-score ceilings passed. The declared recovery milestone was not met. Both component errors remain worse than the original pre-water-fit baseline 13.17/13.10.

The independent actual audit passed 60,601 checks across 218 parameter vectors, including 192 new vectors, 4,968 new native cells and 6,048 programs. It verified the final 46-cell selection. Failed SLSQP statuses at steps 2, 4, 5, 6, 7 remain recorded; actual trial admission and the final audit are separate from optimizer success. The stop was the original eight-step limit. No extension or production promotion follows.

## What the saved three-body test shows

Six new serial native evaluations completed the same eight previously consumed E3/mu3 cycles. An independent raw-input/output and Decimal arithmetic review passed. These are development regressions, not untouched validation.

| Paired Dy minus La three-body energy error | Original baseline | Prior selected fit | Component-aware fit |
|---|---:|---:|---:|
| W8 near | -1.4071 | -5.3119 | -4.9105 |
| W8 far | -1.0393 | -4.6977 | -4.3125 |
| W9 near | -2.2043 | -5.7713 | -5.3514 |
| W9 far | -0.8789 | -4.5501 | -4.1647 |

Signed errors are model minus the saved reference, kcal/mol. All four paired energy and dipole errors improve relative to the prior selected fit, but remain worse than the original baseline. Four of eight individual E3 errors improve relative to the prior fit; seven of eight remain worse than baseline. All eight individual dipole errors worsen slightly relative to the prior fit, while remaining better than baseline. Full cycle values and original reference limitations are in the [independent regression report](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/postselection_regression_v1/independent_actual_review_v1/REPORT.txt).

The terminal paired component errors remain P8=+21.4202/B8=-17.4440 and P9=+24.1285/B9=-20.4574 kcal/mol. Their cancellation makes total contrasts look substantially better than their parts. B contains all nonpair orders; it is not uniquely polarization or charge transfer. Neither these comparisons nor the PBE0/frozen-electronic-state references establish a unique missing mechanism or an absolute binding free energy.

## Optimizer and preservation constraints

One independent conic solution of the saved first-step linearized problem agrees with the successful first SLSQP proposal in q to 1.03e-10. This argues against a material first-step objective error; it does not explain later SLSQP failures or prove global model capacity. Clarabel returned AlmostSolved. The full returned fraction has a small negative epigraph margin; the arithmetic audit passed while strict solver success remained false.

Old energy/force, pair-dipole and susceptibility score ceilings bind at that first-step solution; all 25 individual physical gates have slack. These four score ceilings are development preservation rules, not physical uncertainty limits. Therefore the bounded fit does not isolate lack of functional flexibility from restrictions imposed by that preservation policy.

The original Solved-gated ceiling-removal diagnostic remains unexecuted. A separate prospective saved-only diagnostic is being prepared: construct a feasible epigraph baseline without changing its physical parameter vector, then remove only the four score ceilings in one same-origin local solve. Preserve all physical gates/domain/trust limits and report raw solver status and witness feasibility. This diagnostic will not change the completed fit, release a new fit, open held-out data or admit a model for dynamics.

## Cost and exact evidence

Fit: 7729.982wall seconds, 7770.528controllerCPU+1180.655childuser+1567.813childsystemCPU seconds; 12 single-thread native workers, maximum overlap 12. Audit: 383.424wall/354.052CPU seconds. Serial trimer regression: 15.185wall/14.479totalCPU seconds, six native programs; saved review 3.129wall/2.879CPU seconds. All are truthful authorized login work, zero scheduler allocations; inherited computation is not charged again.

- [Terminal fit](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/run_v1/execution/RESULT.json): 744f1faa.
- [Actual fit review](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/independent_actual_review_v1/REVIEW.json): 0e4851c9; selected artifact fdf42081, parameter key e41fcd3d.
- [Trimer result](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/postselection_regression_v1/run_v1/execution/RESULT.json): d4685a09; [independent review](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/postselection_regression_v1/independent_actual_review_v1/REVIEW.json): 0e519767.
- [Conic arithmetic review](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/conic_saved_comparison_v1/independent_actual_review_v1/actual_v1/REVIEW.json): c59a4c5a, source result 341472cd. No strict solver admission.

Written 22 October 2026, 23:32 UTC. Goal ACTIVE; no qualified interaction model or exchange prediction yet.

## Subsequent local tradeoff result (23:36 UTC)

The separately declared saved-only solve completed once, retaining AlmostSolved status. The unchanged-domain, unchanged-physical-gate full witness is strictly feasible in the linearized arithmetic checks. Its actual normalized component score is **0.438713**, versus **0.842764** for the original full-constraint baseline. Removing the four development ceilings therefore exposes substantial local freedom; it does not establish native accuracy or a global optimum.

The price is explicit: old E/F, pair-dipole, collective-dipole and susceptibility squared-score ratios become **3.0574, 1.0150, 1.0316 and 37.7091** times their old ceilings. All individual physical gates still pass in prediction, with some close to their bounds. This is a tradeoff among observables, not an overall model improvement.

The policy written before this solve selects the first strictly feasible fraction with at least0.05 absolute q gain for one exact nonlinear test. The full fraction qualifies (gain0.40405); independent actual review and its exact one-vector adapter are being prepared. No new fit, second fraction, reference calculation or held-out test is authorized by this conditional policy. Solver cost1.568wall/1.500CPU seconds, one thread, zero native calls; result1d4c5db9.
