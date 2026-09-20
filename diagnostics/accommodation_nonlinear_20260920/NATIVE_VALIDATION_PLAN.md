# Native validation of energy-selected donor candidates

2026-09-20. Declared before nonlinear pilot results. Within Jacob's authorized
overnight accuracy research; this is separate from the first cheap pilot.

## Question and fixed selection

Does native r2SCAN-3c confirm the energy lowering and local forces at the donor
arrangements selected by the cheap composite? This tests the model outside the
previous +/-0.2rad Asp checks, including any coupled Glu response. A cheap minimum
is a minimum in a small frozen-scaffold subspace, not automatically a DFT minimum.

The declared set is all eight endpoint candidates from the first four-case
nonlinear pilot, one final candidate per metal/case. Select only the optimizer's
final candidate; never choose by a class label or a favorable score. A missing
candidate remains unrun, with its reason. Keep available boundary/nonstationary
candidates explicitly marked if others qualify; do not call them minima.

Preparation may proceed once the final immutable pilot collection exists.
Submit only if at least one case has both cheap endpoints pass the pilot's
interior/gradient/curvature criteria. If no pair qualifies, retain a dry-run
manifest and first diagnose that concrete failure. This is a mechanistic gate,
not a requirement for a favorable classification. No default or threshold changes.

## Native calculation and original-state reuse

At most eight new analytic endpoints, using exactly:

    ! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 EnGrad

Preserve original charge/multiplicity, atoms/order, caps, source/mapping and each
final candidate coordinate file. Add EnGrad only to the prior native recipe;
do not silently add TightSCF, a new ECP/basis or numerical gradients. Record
actual SCF tolerances and any automatic EnGrad-related numerical change. Verify
observed native SCF/XC/CPCM, dispersion, gCP and (for La) ECP analytic drivers.

Reuse exact expanded-context q0 energies from the old second_shell manifest
for1H4I/4MAE and the completed original torsion-native manifest for the two PLM
cases. Verify original receipts, states, method and actual executed XYZ; retain
the known sub-femtometre map roundoff separately. No q0 recomputation is needed.
If a compatible origin is absent, the associated work is unavailable, not zero.

Extract the complete Cartesian analytic gradient and project with the physical
Jacobian at the actual candidate q, including cap chain rules. Retain active-mode
derivatives and all inactive-mode restrictions. No Hessian, displacement scan,
native geometry optimization, entropy or nuclear population is requested.

## Comparison and execution

Report actual native work E(q_candidate)-E(q0), MACE/solvent/composite work,
endpoint errors, and the change in Ca-minus-La contrast when both endpoints are
available. Report raw native residual torques beside cheap residual torques.
No new general accuracy tolerance or classification-dependent pass criterion is
introduced. A native energy decrease supports the chosen motion; nonzero native
torque may still disqualify calling the candidate a DFT stationary geometry.
Report physical usefulness, discriminatory usefulness and cost separately.

Use the established ORCA manifest runner,16 ranks per endpoint, up to four
concurrent endpoints per64CPU allocation,600GiB RAM as in the working native
comparison. Split the fixed task list into at most two allocations for throughput.
No project time/core-second budget, outcome-selected reruns or hidden retries.
Preserve every attempted output and count actual new evaluations/allocations.
Existing native torsion, baseline-fold and independent CC jobs remain untouched.

All artifacts belong under workspaces/accommodation_nonlinear_20260920/. The
adapter should only connect existing candidate/receipt/gradient machinery, not
introduce another scientific engine. These results support an opt-in research
decision; biological utility still requires the same rule on known references.
