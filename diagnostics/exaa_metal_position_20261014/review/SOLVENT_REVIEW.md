# Solvent continuation: static implementation clearance

**No launch-blocking issue was found** in frozen `solvent_v1/CONFIG.json` and `SOLVENT.py`. This is source/schema inspection, not a passed molecular evaluation. Root owns Slurm execution and must retain actual preflight/cell failures.

## Energy, reuse and parameters

The adapter reuses all eight completed transplanted MACE values and all eight exact original component records. New endpoints use the unchanged explicit-cell loop in `MATRIX.solvent`, so the former three-candidate `summarize()` is not applied. Each new geometry/metal has standalone vacuum/water evaluations in both xTB and tblite, followed by the existing native water-charge diagonal check and vacuum-charge ALPB replay. The finite maximum is 32 SCC calls plus 16 native component evaluations; no MACE, DFT or optimization call is added.

Original reuse checks exact XYZ/state/MACE/runtime joins, re-evaluates the saved-artifact xTB/tblite **component-v2 parity**, and retains native water-diagonal, cavity-identity and informative charge-rounding checks. The retained absolute-point acceptance flags are reported separately; an old absolute-v1 failure is not relabeled as passing. Original numerical data are read, not recomputed with a different backend.

Executed SCC parameters are the pinned **published standalone GFN2 file**, SHA256 `f6f05c6c5264cb2f75c54f8ac6838044541c1e493b16f1a4d4ee4fa64de5f545`, with pinned xTB executable `481a1e7d8151fc08e0341b8e06ab60ac09cf1c2628738eeb83346271274bf8f4`. Native ALPB uses the already-qualified internal factory/binary. The full corrected ORCA JSON (`a4f83ca72d58d6e80f634e61bb8e99b7af0f26c38cc4a9295e85215ecf72d05c`) is explicitly **audit provenance**, not an executed ORCA input or a globally reordered standalone parameter file. The adapter checks its correction-policy identity; no ORCA GFN2 call is made.

## Comparisons and failure handling

Primary work is computed from each metal/placement's own origin, with MACE and frozen-solvent terms separately converted to kcal/mol. `D = W_Ca − W_La`, per-host transferred-minus-original `delta_D`, and target-minus-control differences retain the declared signs. A missing new cell makes its dependent response unavailable; the vacuum result is not returned as a successful solvent fallback.

Secondary mathematical minima use the same compatible pool for both metals. The explicit tie order is original origin, original adaptive-La, transferred origin, transferred adaptive-La. The original two-geometry pool and expanded four-geometry pool require all their declared metal/geometry cells. Selected XYZ records, per-metal works relative to the original origin and the change in differential work remain visible. There is no operational 0.1 kcal rule, classifier band, new calibration or population interpretation. This pool lacks historical adaptive-Ca and is not the full historical classifier pool.

The secondary analysis was declared after the vacuum result and before solvent execution, as documented in `SOLVENT_PLAN.md`; it remains development evidence. A complete mathematical pool does not demonstrate a physically optimized metal position or improved discrimination.

## Frozen records

- Config: `workspaces/exaa_metal_position_20261014/solvent_v1/CONFIG.json`, SHA256 `0f162cb433686011436a4c9f04a71ebded1e6c3e8de112383032e7ec6933d708`.
- Adapter: `cc02ca883ec2f8504e007674b2f6273d4413a2361d0d1f555efe00888bf9a3c8`.
- Plan: `790eb7f4ec0af3b163cf4b9b4ce6c9d65f145644887a2954f3afb53c69f4054d`.
- Native manifest: `bcd9fc9d8923139dfcb8d7b8d0191c5c29484e700ddc8ffd488c1dea09c10b41`; replay binary: `3d6a9d23aa1d1152c368b21f7f6d55a135279596b8d01a3a1b3ebc42617e75b6`.

The real Slurm guard and seven-CPU allocation precede scientific imports/execution. Existing endpoints use at most seven independent one-thread workers; libraries/build/runtime are pinned. This reviewer performed no numerical scientific execution, tests, submissions, threshold changes or edits to the frozen implementation.
