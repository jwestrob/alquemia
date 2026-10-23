# Fresh CPU-refactor numerical qualification — passed

The [independent saved-output review](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/sampling_cpu_refactor_v1/independent_actual_review_v1/actual_v1/REVIEW.json) passed all57 fresh programs,
636 internal finite-difference energies and38,992 recorded checks. All four
required qualification fields are true. The result is bound to source9ccee436,
buildf1bee775, runner90616207 and actual result5b3a21f2; executable hashes are
retained in the review. No old27-output prefix was reused.

- All21 original energy/force identity comparisons are exact.
- Maximum Cartesian derivative discrepancy:8e−8 kcal/mol/Å, below the original1e−3 gate.
- Maximum rigid radial/rotation derivative discrepancy:6.88166e−5 kcal/mol/Å; both step sizes pass unchanged printing/derivative bounds.
- All12 raw matrix snapshots, both RHS and dual solves, reconstructed H/h, stationary polarization energies and zero-field energy-route checks pass. Analyze/testgrad energy-route discrepancy is zero.
- All703 induction summaries are admitted; maximum printed residual9.9e−9 Debye.
- Normalized spectral endpoints across snapshots span0.611100629–1.815593274; each passes the original relative SPD criterion.
- Maximum odd-field conjugate-moment discrepancy:4.33447e−9 ea0; even-field/matrix-curvature discrepancy:1.06075e−5 a0³. Both pass original condition/printing bounds.

All raw values replay exactly. No recomputed arithmetic difference was observed;
the prospectively declared derived-value allowances were therefore unused.

Actual serial qualification cost:95.383s wall, 91.063021 measured controller+child
CPU-seconds, zero scheduler allocation,1 native worker plus1 controller. This
saved review cost8.762s wall/8.053CPU-seconds on1 numerical thread, with zero
native calls or reruns. Historical/inherited costs are not added.

This admits numerical consistency only for the declared vacuum single-metal
water fixtures and fixed profiles. It does not prove physical accuracy, global
stability, periodic/GPU support or production sampling. Old/new functional parity
is a separate review; no fitting or held-out outcome was inspected here.
