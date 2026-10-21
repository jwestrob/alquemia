# Independent actual preparation review:16667

21 October 2026 UTC. **PASS: actual Slurm16667 completed the declared preparation with status `prepared` and zero engine calls.** No input/preparation blocker was found to the separately released resource-only endpoint. This review used saved metadata, hashes and literal byte/field equality; no scientific replay or numerical imports were performed on login.

## Exact evidence

| Artifact | SHA256 |
| --- | --- |
| preparation/execution/RESULT.json | `cd2d06b222beceb15eaa03a32e5802d41a7334215241b5bf45cf4be7948484f5` |
| INPUT_EVIDENCE.json | `67bc9bca6dd7ffd7e23d667a81dfeda83238f23355e44020bc4a17cefb365d11` |
| MANIFEST.json | `d0c662efac57ed00ee196fac1502ddd07e8487d3c9f467f20cc89e2ebd3e2082` |
| TARGET.py | `5cad26238cc71e676a5af341ea8e0fe2a3ef4d7476f3a3d0503879c27ecc35c3` |
| CONFIG_PREPARE.json | `32c8501b73d00df7affa6bf0a6d86937fe0ca2cf0b0028f36c5aa6bb37f96c9a` |
| central_La/endpoint.inp | `4359537c4a569f2a68793f60a509bda509c854a15a7aaf1060790a9407086896` |
| central_La/physical.xyz | `55e07d8146adc13c09203243849b65e80a6afef89b181cfe4c22b265a994b72d` |
| central_La/def2-tzvpd.bas | `7389c56f0cff9e76098b0757a34d1abdd79c48e9249597cea074a0106a3462a6` |

All selected pins match. Input, XYZ and external-basis bytes are exactly identical to16611. Its physical-parent, positions, complete atom-ID ordering, atom count, target index, charge/multiplicities and electronic-state fields are unchanged. The sole manifest task equals TASK.json; configuration/result/manifest scientific bindings agree. Input evidence preserves all nine OOM16611/original16602 parent pins without substitution.

The actual on-Slurm preparation passed the frozen source/physical/input gates. Saved state remains whole256 La at index255, q−2 UKS singlet,1030 explicit electrons,515/515 occupations, external-file ECP46, PModel and beta rotation `{514,515,20,1,1}`; expected orbital dimension remains6239. The prepared scientific deck is unchanged, including the EnGrad request that selects SCNL in this version.

The new manifest/configuration consistently declare **12MPI ranks,24 allocated CPUs,two CPUs per task and one worker**. This is the future run plan, not an allocation observed during the one-CPU preparation. Actual96000MiB, single-thread execution and derived MaxCore6200 must still pass the run's resource guard. Preparation does not prove adequate memory or convergence.

All runner/DFT/MACE/ALPB/GFN/MM, geometry, optimizer, displaced and other-metal counts are zero. Native `unavailable/not_attempted`, comparison unavailable and null atomic-guess artifacts are appropriate for preparation. The task directory contains only TASK.json, input, XYZ and basis. Application elapsed time is5.926165568875149s; actual scheduler cost remains separate.

Root must bind this actual result and its single-endpoint decision before execution. Retain both failures, fresh unseeded files, native atomic-initialization accounting, unchanged scientific admission and pending method/state/gradient review. No automatic escalation, solver change, displaced pair, optimization or accuracy/physical-state claim follows from this pass.
