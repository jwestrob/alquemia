# Constrained proposal: optimizer failure, saved feasible half step

The single declared SLSQP solve returned status 8, “Positive directional
derivative for linesearch,” after 77 iterations. **The original proposal run is
failed under its unchanged admission rule**, which requires solver success.
There is no repeat solve, tolerance relaxation or claim of convergence.

The result retains all predictions at zero, half and full returned step. Both
nonzero points pass the 25 original predicted energy/force gates exactly. The
full point has a maximum normalized constraint excess of 3.48e-9 and is not
strictly feasible. The predeclared half point is strictly feasible in the stored
linearized problem: its smallest inequality margin is +1.297e-6.

| RMS component error, kcal/mol | Origin | Half step | Full step |
| --- | ---: | ---: | ---: |
| Pair contribution | 26.6824 | 25.5887 | 24.4950 |
| Collective remainder | 22.3032 | 21.3889 | 20.4748 |

All four old block norms decrease at half step. These are predictions from a
local model of the residuals, not native energies or forces. The large remaining
errors and all reference limitations remain.

The independent arithmetic check passed all 851 inequalities with strictly
positive margins and all 25 original gates; scalar replay differences are at
most 3.6e-15. Its receipt is the matching workspace
`independent_half_review_v1/REVIEW.json` (194344cb). It establishes feasibility
only in the saved linear model. The original failed optimization stays failed.
Root has separately declared one fixed-half native development pilot in
`../nonlinear_check_v1/ROOT_DECLARATION.json`; source preparation/review and a
separate execution release remain required. This is consumed development, not
a repaired model, optimizer convergence or fresh validation.

The source/preparation passed 12 checks, including 25 physical-residual
reconstructions and 851 analytic constraint derivatives. Independent source
review `a750aeef` passed. The actual solve cost 14.648 seconds wall and 14.600
CPU-seconds on one thread; preparation cost 3.598 wall/3.576 CPU-seconds.
There were zero native/reference calls or scheduler allocations.

[Actual result](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/constrained_step_v1/actual_v1/RESULT.json)
and [unchanged specification](SPEC.md) retain the failure and exact constraints.
