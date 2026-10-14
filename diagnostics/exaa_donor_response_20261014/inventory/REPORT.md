# Exact ExaA donor-response inputs

The original **192-atom** contexts are available for the two-Ca sample3 target and sample0 control, with complete saved origin, Ca-proposed and La-proposed geometry pools for both metals. No exactly compatible completed native DFT endpoint or wavefunction was identified in the named prior collections. This inventory supplies inputs for the next reference test; it contains no new molecular result.

## Runnable input inventory

- [INVENTORY.json](../../../workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json), SHA256 `757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343`.
- `rows[0]`: `q9z4j7-pqq-la_model__conditioned_2Ca_v1__seed-1_sample-3`, target.
- `rows[1]`: `q9z4j7-pqq-la_model__conditioned_2Ca_v1__seed-1_sample-0`, control.
- Each row exposes `matrix[Ca|La][origin|adaptive_Ca|adaptive_La]`: exact XYZ, charge, multiplicity, atom count, saved MACE receipt/request/model, native frozen-ALPB receipt and components, underlying GFN2 inputs, and previously computed origin-relative works. All **12 cells** remain visible, including the Ca proposals outside the initial eight-cell reference plan.
- Maps, source preparation, full-source and spectator admission, state signature and original optimizer result pins are retained per row. `BUILD.py` performs metadata/hash joins only; numerical verification belongs in the Slurm reference preparation.

The shared state signature is `63469aa8c373d6c2d3951e2f69ff7c51dfb9bf108a28b1921800a23078393532`: Ca charge −3, La charge −2, physical multiplicity1, dry oxidized-PQQ3− preparation hypothesis, no explicit coordinating waters. The archived GFN2 valence-electron count588 is **not** a DFT all-electron/ECP electron-count prescription. Both native Ca atoms remain in the full source; the spectator lies outside the local context. Target sample3 uses chainD as the PQQ-associated Ca, while control sample0 uses chainB.

These are the original contexts preserved by the [packing manifest](../../../workspaces/exaa_source_packing_20261014/packing_v1/manifest.json), not the subsequently expanded 249-atom packing matrices. No geometry, proton inventory, cap or state was changed for this inventory.

## Saved response to be checked

Values below are copied from the already executed [small-context comparison](../../../workspaces/exaa_source_packing_20261014/packing_comparison_v1/execution/RESULT.json). Work means candidate energy minus that metal's origin energy, in kcal/mol. The current candidate is MACE vacuum energy plus native ALPB evaluated on the GFN2 vacuum charge distribution; it is not the older self-consistent GFN2 solvent subtraction.

| Source | Ca work at Ca proposal | Ca work at La proposal | La work at Ca proposal | La work at La proposal |
|---|---:|---:|---:|---:|
| Target sample3 | −9.698819 | −7.898852 | −25.271250 | −26.390560 |
| Control sample0 | −7.624675 | −6.395494 | −12.326184 | −12.934776 |

Both operational and mathematical selection choose the Ca proposal for Ca and the La proposal for La. The selected Ca-minus-La differential work is +16.691741 for sample3 and +5.310100 for sample0. The target's selected differential decomposes into +15.495668 MACE and +1.196073 frozen-ALPB; the control gives +4.650659 and +0.659441. Full unrounded quantities and all components are preserved in the inventory. These are local model responses, not affinities or new class-validation results.

The parent-owned first reference stage uses origin and `adaptive_La` for each metal in both sources, retaining the Ca-proposed coordinates for a possible subsequent complete-pool comparison. The inventory does not prescribe a new electronic Hamiltonian or import calibration bands into DFT response checks.

## Existing DFT and wavefunctions

The pinned [preserved DFT collection](../../../workspaces/accommodation_fold_DFT_20260920/run_v1/final_collection.json), SHA256 `9d7ae05936a28cd845ae63f2da34384218f1f23f3968da7909a346235e1ec559`, contains ten historical Q9 source rows:

- Five original single-Ca folds were unsupported at preparation (`exact selected PQQ required`). They provide no completed endpoint.
- Five La-conditioned folds have both Ca and La origin references: **73 atoms**, charges −2/−1, singlets, r2SCAN-3c/CPCM(Water)/DefGrid3. Actual `.gbw` files remain beside these outputs and are pinned as historical directory contents. The canonical directory contains both metals' wavefunctions; directory membership is not a current-endpoint restart association.

Those references differ in source conditioning, atom membership, charge and geometry from the current two-Ca 192-atom matrices. Neither their energies nor their wavefunctions are exact reusable endpoints here. The [earlier bounded Q9 audit](../../orca_parameter_audit_20261012/q9_response_diagnosis/DFT_INVENTORY.json) also identified no matching donor-motion reference. DFT checks on other proteins' donor motions do not qualify these ExaA proposals.

This is a search of the named authoritative collection, prior donor-response audit and current two-Ca source trees, not a claim to have searched every external archive. Native ALPB replay output is not native DFT. Exact-compatible DFT cells and wavefunctions are explicitly empty in the inventory.

## Execution status

Metadata preparation completed; zero molecular calls, zero scientific tests and zero submissions. Root owns the new reference plan and its Slurm execution. The existing production scorer and all archived evidence remain unchanged.
