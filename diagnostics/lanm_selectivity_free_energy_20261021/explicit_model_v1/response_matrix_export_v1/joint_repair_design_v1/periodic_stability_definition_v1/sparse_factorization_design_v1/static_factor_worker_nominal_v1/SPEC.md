# One nominal static factor attempt — source only

This version applies NOMINAL_CONTRACT `2aa3ac43` to the preserved worker `04624424`. It changes only nominal-shift arithmetic, explicit scope/status, and input provenance. No assembly, factorization, native call, or solve is authorized by preparation. The parent worker and all its outputs remain unchanged.

For the one saved La BIO operator, assemble the unshifted normalized real+self+candidate stiffness matrix `S=B0+H`. Use `U_nominal=U_B+U_R` and `tau=1e-6*U_nominal+delta_odd`; factor `S-tau I`. `U_R` already uses raw maximum absolute q, so the odd-grid term is added to the lower-margin shift only. `e_app` remains **null**, with `application_error_bound_established:false`; null is never interpreted as an established zero error. The original relative threshold is unchanged.

The new context joins the exact composition, terminal candidate/review, initialized La PME export, full-box execution, independent review, and root nominal-data admission. `U_R=4900.1575968087745` and `delta_odd=2.3547987346895603e-12` must match all three saved reports from that snapshot. Their admission is ordinary-double numerical evidence; the export's `PSD_admitted:false` is preserved. The review flag is consequently `initialized_PME_nominal_data_and_bound_reviewed`, not a claim of rigorous PSD/application proof.

Matrix arrays, `U_B`, producer/result pins, root input review, source review of the final bound configuration, and root release remain unbound. The existing 35874-dimensional La graph and 16371009 lower entries remain exact; `.npy` data are float64 and both index arrays int32. All sites have positive intrinsic alpha and all mutual scales are one. A lower triangle alone does not establish full-native symmetry; the upstream review still owns that check under the unchanged 1e-10 allowance.

The qualified backend and factor path are unchanged: one symbolic analysis, one supernodal numerical attempt, explicit converted LL storage, and at most ten sampled sparse action/solve checks. No alternate ordering, retry, fallback, rank update, or additional probe is introduced. The unchanged isolated child envelope is 32 GiB address space, 900 s wall, 10800 CPU-s, OpenBLAS 12/OpenMP 1, and at least 64 GiB host MemAvailable before start. Host headroom is not a reservation. No factor-array import/read occurs before child limits.

Success is `passed_nominal_static_sufficient_numerical_test`, limited to this saved operator and the existing ordinary-double sampled checks. Failure of the sufficient real split remains inconclusive for full PME stability; malformed inputs, backend errors, warnings, and time/memory limits retain their existing separate classifications. Neither outcome admits trajectories, periodic response forces, physical accuracy, or MD.

## Finite source checks declared before execution

Run TEST_SOURCE.py once, stdlib only, one thread: nominal formula and unknown-e_app behavior; finite/sign and mismatched-value rejection; synthetic exact PME/candidate/export/geometry/producer joins and deliberate mismatches; exact supervisor, resource, factor, and audit call preservation. Compile the source without execution. No backend imports, matrix reads, geometric calculations, or molecular calls. Save all checks and wall/CPU cost even on failure.
