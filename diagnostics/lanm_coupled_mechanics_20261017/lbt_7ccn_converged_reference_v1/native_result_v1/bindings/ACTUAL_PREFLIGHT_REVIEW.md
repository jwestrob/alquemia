# Independent actual preflight review — 15543

**PASS; no blocker to binding the declared single native run.** Actual status is `preflight_passed`, with frozen source `f251ca1e` and CONFIG_PREFLIGHT `346f9622`. Review used saved metadata, exact JSON/deck equality and 26 selected file-pin checks; no scientific import or numerical replay occurred.

- The manifest contains exactly one `new_La` analytic-gradient task: 256 physical atoms, target index255, charge −2, physical/effective singlet, ECP46 and 1,030 explicit electrons. Execution resources remain 12 MPI ranks and one worker. Task state, physical parent, positions, physical IDs and cached ALPB match the exact admitted step725 terminal POINT `6605f43c`.
- The native input is byte-identical to the actual step92 reference template `b94c22fc`. No new native output or gradient file exists at preflight review.
- The reused `old_La` row equals actual 15321's original `new_La` row and collector 15456's corresponding row exactly, including native task identity, receipt and residual observations. Both cached solvents are available. The original source/recipe, step725 checkpoint/search and physical guard match job15525; source rebasing remains false.
- All molecular-call, native-runner, optimizer-step, new-geometry and other-metal counters are zero. Selected tasks remain empty. Only the old endpoint has an available native row/force report; the new row remains unavailable, `full_two_available=false`, and comparison is unavailable. Preflight does not admit a future reference result.
- The preflight reports step725-minus92 model work of MACE **+16.362172**, fixed ALPB **−45.442067**, composite **−29.079895 kcal/mol**. These are saved producer results; this review did not recalculate them or substitute the shorter job15525 continuation work.

The old reference retains EnergyCheck admission: energy change passes its printed tolerance while density, orbital-gradient, orbital-rotation and DIIS flags fail. Physical-state and force-consistency qualification remain false. This preflight does not establish reference stationarity, physical accuracy or solvent validity.

Root must bind this exact actual preflight and the execution decision before the sole new 12-MPI call. Preserve the frozen scientific binding, native settings/admission, unavailable outcomes and no-retry scope.

## Exact pins

- RESULT: `5edd76c0b7826b21eb2d92f79f7989355322496f351aa04cde8d70fa0554fc82`.
- MANIFEST: `e1e786d3e2b5c2284fad2248378dd81413049a8b281a6521d7ed6fa656dbb2ee`.
- DIRECTIONS: `2339b196a609e98c6c27b4c70e3d86f3ffe212bd6ef72cc2f8951d482bce7305`.
- Native input: `b94c22fc46ff47bd8c678f2e49634948fa99195638f6f69823d0b99f075b102d`.
- Prepared XYZ: `780454c2284989ed2f2c27633102056168bd33e27c22f3c34bc859525ab576b4`.
- Reused native receipt: `6583718c6fbe8aa5aa1843b2ab1fefc9e877977c9990f484560f00a04be1b9d9`.

Actual evidence is under `workspaces/lanm_coupled_mechanics_20261017/lbt_7ccn_converged_reference_v1/preflight/execution/`. Application elapsed/allocated-CPU time is 5.567950 seconds; scheduler accounting remains separate.
