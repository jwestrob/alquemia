# Matching La/Dy large-core reference assets

**Matching La orbital/ECP and a real exported La AutoAuxJ asset are available.** Job1220299 completed on standard-shared with oneCPU/mem0, four allocatedCPU-seconds. This was basis construction only: no molecular SCF iterations, atomic fitting SCF, final energy, gradient, or selectivity result. Production is unchanged.

## Orbital/ECP provenance

| Asset | SHA256 |
|---|---|
| `legacy/qmmm_lc/ecp_lib/lcecp1_tzvp_La.bse` | `4aa9b5c40669bb8b17f0641b1df65e333edd7714dd89154d0ea8f38fa304d266` |
| `legacy/qmmm_lc/ecp_lib/orca_La.ecp_basis` | `c5298b309b2d897e5de59ddb6b4c7992913dd8cba18b7cd21f9757ce09b167d3` |
| `benchmarks/hans_lanm_dy_qmmm_correction_v1/qmmm_calibration/basis/lcecp1_tzvp_Dy.orca.inc` | `87bb209a42eca365c93d8f69088a57371285623a1d4afa1c1336390279fafb20` |
| Same directory, `lcecp1_tzvp_Dy_autoauxj_orca611.inc` | `ccefb184b4504f1dc9fb1a2c8d49656ff5e0de4bc5637c18db9b5d495900b8d3` |
| This directory, `lcecp1_tzvp_La_autoauxj_orca611.inc` | `31d3271de21c3c961f8f6c48b9184ab43c57db8aa0c9c6c58a259e5741301d81` |

The La BSE file identifies lcecp-1-TZVP version1, data from authors, BSE0.12. Reapplying the existing converter reproduces `orca_La.ecp_basis` **byte exactly**. Conversion rounds orbital numbers to ten significant digits: maximum absolute/relative numerical difference from BSE is4.999997e−9/4.526627e−10. This rounding is disclosed, not asserted to be bit-identical to the original numeric source. ECP numeric text is retained. The source is the occupation-specific family in the [primary basis paper](https://pubs.rsc.org/en/content/articlehtml/2026/cp/d5cp04944j), not the native ORCA def2-La basis merely because both remove46electrons.

La orbital contraction is9s8p5d2f→6s4p3d2f,47spherical functions. The real export confirms ECP46 and these functions. Keep the original orbital coefficients: printed ORCA orbital coefficients are normalized and must not silently replace the pinned input asset.

## States and electron bookkeeping

| Element | Physical trivalent hypothesis | Frozen core | Explicit electrons, neutral/trivalent | Valence backend |
|---|---|---:|---:|---|
| La | 4f0, physical singlet |46|11/8|restricted singlet|
| Dy | 4f9, physical sextet |55|11/8|restricted effective singlet|

At identical ligand atoms and charge−1, the consumed50-atom core therefore has208explicit electrons for **both** La and Dy (104alpha/104beta). The previously prepared195-atom Hans region would have806explicit electrons in both large-core models if reused unchanged; that is bookkeeping, not a completed calculation. The native explicit-f Dy count833 belongs to a different representation.

For any new region recompute sum(Z−Ncore)−charge from its actual atoms. The physical Dy sextet metadata stays distinct from its effective closed-shell valence solver convention. No artificial change of physical spin is being made to make the original model converge.

## Actual auxiliary export

Workspace: `workspaces/metal_environment_response_20260926/la_aux_export_v1/`. `export.inp`, runtime input, `export.out`, `receipt.json`, accounting, allocation module and submission script are preserved and pinned in LADY_REFERENCE_ASSETS.json.

ORCA6.1.1 ran `PBE0 RKS RIJCOSX def2/J PrintBasis NoAutoStart`, explicit `%scf DryRun true; Guess HCore; AutoStart false` (as multiline block), the pinned La orbital/ECP, then element-specific `NewAuxJGTO La "AutoAux" end` with AutoAuxSize1. The atomic coordinate origin and La(III) charge+3/singlet are **basis-construction input**, not fabricated molecular evidence. The memory setting was derived inside the allocation from the existing resource policy.

The actual output stops after the SCF memory estimate, with normal termination. It contains HCore initialization but no SCF iteration table, no atomic density fitting subprocess, and no final single-point energy. ORCA wall time was3.23920655s; Slurm elapsed4s/1CPU/zeroGPU. This differs from the old Dy NoIter/PModel export, which did invoke an internal atomic fitting SCF. No fallback was needed.

Only the printed `NewAuxJGTO La` block was extracted. Its65uncontracted shells are15s12p11d10f9g5h3i,351spherical functions; every exponent is positive and coefficient1. The terminator was normalized from `end;` to `end`. Exponents/coefficients retain the exact printed values. The literal asset has no enclosing `%basis` wrapper. No energy zero from a dry run was recorded as science.

## Use and remaining qualification

For the declared conventional PBE0-D4 reference, retain def2-TZVP and def2/J on light atoms; install each metal's own matched lcecp-1-TZVP/ECP and literal AutoAuxJ inside one `%basis` block. This follows the [ORCA6.1 AutoAux interface](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basisset.html). Do not add generic AutoAux at runtime, paste Dy exponents onto La, substitute native La basis, or mix native-3c reference offsets. The export proves basis-generation capability, not fitting-energy accuracy.

A matched La force endpoint and balanced La/Dy response comparison remain new scientific work owned by root. The ongoing Dy force scout and this export do not establish affinity, basis convergence, auxiliary accuracy or full hybrid force qualification. No reference score or production setting was changed.
