# Complete finite gradient assembly demonstrated on real matching fixtures

**Both archived La A/B electronic gradients now combine with actual matching classical gradients on all1,891physical atoms.** The two authorized classical evaluations completed in job1220328: one allocatedCPU-second, zeroGPU, no new quantum calculations. This closes the missing assembly operation, not the scientific qualification of a full protein response model.

## What ran and what was matched

Eight OpenMM Reference energy/force queries were executed: exactly two original physical geometries × four components. Before any query, preparation verified that the selected repaired-source ledger and the old fixtures have identical native System hash, full physical inventory, QM membership, shifted point charges and cap definitions. Each geometry was replaced explicitly with its archived physical A/B coordinates. No topology, charge, parameter, region, cap or boundary choice changed with that substitution.

The selected ledger names the target `target_EF3`; the assembly names its real source site `A/203//LA`. An explicit one-to-one alias maps them after matching source identity/position. This is not a second atom or a changed metal. Other source IDs match directly; the assembly reorders classical rows by ID.

The evaluated terms were retained native bonded, MM–MM LJ, shifted-charge MM–MM Coulomb and realQM–MM LJ, using the previously selected boundary exclusions and actual pure12–6 ion parameters. No capFF particle, classicalQM–MM Coulomb, QM–QM classical term or C4 induction was added. Each component's actual Cartesian gradient and energy was written **before** summation.

| Component energy, kcal/mol | Original A | Original B |
|---|---:|---:|
| Retained bonded |5374.829818064|5374.829818064|
| MM–MM LJ |52.355168811|52.355168811|
| MM–MM Coulomb |−4016.601581431|−4016.601581431|
| RealQM–MM LJ |−56.039607809|−56.037607406|

These are components of the declared finite model, not affinities. The electronic contribution is the archived **native r2SCAN-3c original-H** endpoint, not the repaired frozen-f endpoint. No old energy fills any unavailable new scout.

## Actual force accounting

The electronic component includes realQM atom gradients, J_retained^T times each cap gradient, J_omitted^T times each cap gradient, and all serialized point-charge gradients. The classical component is already a physical Cartesian gradient. Both use kcal/mol/angstrom after one native-unit conversion. The total is their sum; force is its negative.

Actual total translation-gradient residual norms are5.53479e−6 and5.57199e−6kcal/mol/angstrom for A/B, consistent with the saved electronic residual. No force correction was subtracted to force this sum to zero. Both cap-anchor paths are tested, including uncharged omittedCA anchors.

Six real-artifact tests pass. They include the actual classical addition under different row order, exact sign, and explicit rejection of a deliberately corrupted copy of a real coordinate record. Existing mapping tests exercise real saved cap gradients and geometry, missing endpoints and mismatched source states. No fabricated scientific energies or forces were substituted.

## Artifacts and next use

- `workspaces/metal_environment_response_20260926/force_assembly_classical_v1/manifest.json`: exact preparation/execution pins.
- `{A,B}/{retained_bonded,MM_LJ,MM_Coulomb,QM_MM_LJ}.json`: actual per-component energies and Cartesian gradients.
- `{A,B}/CLASSICAL.json`: compatible classical bridge with IDs, coordinates and exclusions.
- `{A,B}/ASSEMBLED.json`: actual mapped electronic, classical and total physical gradients/forces.
- `RECEIPT.json`, `ACCOUNTING.txt`: actual job/cost evidence.
- `COMPLETE_ASSEMBLY_RESULT.json`: compact joinable summary and artifact hashes.

The reusable CLI is `scripts/metal_environment_force_assembly.py` with matching `--inputs`, `--configuration`, `--metal`, `--collection`, `--task-id`, optional `--classical` and fresh `--output`. Supplying classical data triggers exact physical-ID/coordinate compatibility checks; omitting it leaves total gradient/force null. The completed two-cell classical manifest must not be rerun as an automatic continuation.

## Limits remain explicit

The earlier coarse and half-step classical finite-difference gates remain failed; their diagnosed LJ truncation behavior is not a retroactive pass. This assembly test adds no independent full electronic+classical directional derivative check. Fixed-cap double-counting conventions, unqualified ion cross parameters, electronic numerical/state uncertainty and missing solvent remain model limitations. The finite dry sum is executable, but it does not qualify a solvent-consistent whole-protein relaxation, affinity or global LanM preference calculation. A future repaired frozen-f gradient requires its own exactly matching input/field/classical geometry; the module will not silently mix those states.
