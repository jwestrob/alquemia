# Independent actual preparation review

21 October 2026 UTC. **PASS: actual Slurm16601 is `prepared`; no preparation blocker to the separately released, single24-MPI central EnGrad.** Review was limited to saved JSON/text, literal identity checks and hashes. No scientific imports, vector arithmetic, molecular calls or replay were performed on login.

## Exact evidence

| Artifact | SHA256 |
| --- | --- |
| preparation/execution/RESULT.json | `6c4200a5213a820c8acda4f20f62b1f8c4377c2205738a57f41ef182b07117da` |
| INPUT_EVIDENCE.json | `a963fc0fe57b398e2c68da129b6978a86afe0b2845f931f068c2db4d8264ceeb` |
| MANIFEST.json | `84f4b36b1c58fe130abdd90b35989e5a8225ca92fa813e0591d893a55e2c86c1` |
| TARGET.py | `aa62db4f37de3c50765ad81115dd75cb4aef2a99d92b4465d9bbe59df9a4a74f` |
| CONFIG_PREPARE.json | `f8a105695715b69cd0e8fe06ff10b5a032c55d0a97ace77745bc3b6c75352e74` |
| central_La/endpoint.inp | `d482826aa2384c0e00f8c8fb0cfa5dc965e22a3c2a99f17a4e608602a45200c6` |
| central_La/physical.xyz | `55e07d8146adc13c09203243849b65e80a6afef89b181cfe4c22b265a994b72d` |
| central_La/def2-tzvpd.bas | `7389c56f0cff9e76098b0757a34d1abdd79c48e9249597cea074a0106a3462a6` |

The selected artifact hashes, including the pinned ORCA executable, match. Configuration, result and manifest scientific bindings agree. The sole manifest task exactly matches its saved TASK.json; physical-parent, positions, all256 unique physical IDs and electronic-state records agree with INPUT_EVIDENCE. Its old central authority is15545 result `a53125b9…`; the saved force authority is `e6c0d039…`.

The prepared deck retains the published ωB97M-V/def2-TZVPD controls, external GTOName, DIIS/NOSOSCF, Tight/300, DEFGRID3, thresholds, UKS/PAtom/NoAutostart and beta rotation `{514,515,20,1,1}`. The adjacent basis is byte-identical to the frozen public asset. Saved preparation reports6239 orbital functions, charge−2, La at index255, ECP46,1030 explicit electrons and515 alpha/515 beta; actual executed assignments remain to be inspected. The manifest declares one worker and24 MPI ranks.

All engine/runner, optimizer, displaced and other-metal counts are zero. The task directory contains only TASK.json, input, XYZ and basis. The `unavailable/not_attempted` native row and unavailable comparison are appropriate for preparation, not failed chemistry. Application elapsed time is5.965493321185932s; root reports scheduler7s/1CPU, separately from native costs.

This passes saved input preparation only. Native core admission, actual functional/basis/auxiliary/ECP/rotation/VV10 evidence, residual/state diagnostics and measured cost still require review. ORCA6.1.1 is a version-limited recipe port; physical-state, force-accuracy and force-consistency qualifications remain false. No displaced release, retry, optimization or change to the failed PBE0 gate follows from this result.
