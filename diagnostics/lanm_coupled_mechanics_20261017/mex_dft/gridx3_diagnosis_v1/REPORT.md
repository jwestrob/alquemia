# Mex force mismatch: saved final-integration diagnosis

17 October 2026 UTC. Saved-data job **13277 completed**, 2 scheduler seconds on one CPU, zero molecular calls. The original failed qualification remains failed.

## Finding

The geometry-dependent final COSX exchange correction accounts for most of the observed energy/gradient mismatch. Direct replay of ORCA's raw point-charge gradient reproduces our assembled physical gradient and directional projection exactly. This argues against a mapping, sign or unit-conversion defect in this tested component.

All entries below are directional energy derivatives, in kcal/mol/Å, along the original external-charge displacement of A/97//O. The step remains ±0.001 Å.

| Saved quantity | Value |
| --- | ---: |
| Central native analytic derivative | -3.673753903 |
| Original final-energy finite difference | -4.087955022 |
| Original finite-difference minus analytic residual | -0.414201119 |
| Last SCF-cycle energy finite difference | -3.684878152 |
| Printed final COSX correction contribution | -0.401606063 |
| Net post-SCF energy-shift contribution | -0.403077118 |
| Residual after subtracting printed COSX contribution | -0.012595056 |
| Residual after subtracting net post-SCF shift | -0.011124001 |

Both diagnostic residuals remain larger than the original 0.005409 tolerance. Neither diagnostic subtraction is an admitted replacement energy definition. The original PBE0 and total derivative gates still fail; MM and fixed-charge ALPB still pass.

## Precision and interpretation

The net post-SCF shift is the printed post-grid total energy minus the last SCF-cycle energy. It must not be equated entirely with COSX: their directional difference is -0.001471054 kcal/mol/Å, beyond the printed rounding bound of 0.000313758. Central and minus scalar closures also exceed their rounding bounds. These unexplained remainders are retained in RESULT.json.

The fixed-core D4 values are identical at printed precision. The final-energy-minus-post-grid derivative is 0.000000248 kcal/mol/Å, within its printed rounding bound. Printed rounding bounds do not bound SCF error. The central/plus/minus runs have different stopping routes and residual flags; the central calculation also used a different historical MPI allocation. No step-halving, electronic-stability or global-force qualification follows from this diagnosis.

ORCA documents using a final COSX grid for final energies and gradients, and includes grid-weight/overlap-fitting derivatives in its analytic implementation. The documented general implementation does not establish how every post-SCF density-response contribution enters external-charge gradients. A density converged on a different grid could therefore explain this pattern, but that is an inference, not a demonstrated implementation defect. See [official COSX documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/RI.html#cosx-numerical-grids) and the separate DOCUMENTATION_REVIEW.json with supporting sources and limitations.

## Decision

Preserve collector13271 and its failed gate. The original remaining six qualification points and dependent eight-cell response panel stay stopped. A separately specified, bounded consistent-grid test is being considered; no new molecular configuration or submission follows from this saved diagnosis. Whole-peptide reference13051 continues independently.

An informative next qualification would compare a matching central gradient and displaced energies with the same SCF and final integration grids, rather than treating a post-hoc subtraction as a correction. The exact protocol, native controls, cost and stopping rules must be recorded before execution.

## Reproduction and evidence

DIAGNOSE.py, SPEC.md and STATIC_REVIEW.json are the reviewed source contract. The result archive contains exact CONFIG.json, DECLARATION.json, SUBMISSION.json, launcher, RESULT.json, scheduler receipt and logs. SUBMISSION.json records the complete executed command and completion watcher. ROOT_REVIEW.json pins the actual result and original failed collector. All original raw native outputs remain pinned through those records.

The saved replay exactly reproduces the original final-energy, MM, ALPB and total finite differences. Raw point-charge row832 maps once to physical atom1019; no additional charge factor or sign reversal is applied. All six unrun endpoint slots remain explicit. Diagnostic execution success establishes only successful analysis of the existing failure.
