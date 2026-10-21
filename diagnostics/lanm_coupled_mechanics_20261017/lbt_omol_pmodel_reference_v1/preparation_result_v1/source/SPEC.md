# Single PModel recovery of the central published-recipe port

21 October 2026 UTC. Root authorized one separately versioned recovery after job16602 failed during PAtom initialization for La. Original source, declaration, input and failed result remain immutable. Source basis: [guess review](../lbt_omol_target_reference_v1/GUESS_RECOVERY_SOURCE_REVIEW.md) and [actual SCNL discovery](../lbt_omol_target_reference_v1/NATIVE_NL_METHOD_DISCOVERY.md).

## Question and sole input change

Can the unchanged step725 whole256 La endpoint complete with PModel initialization, permitting the declared central MACE/reference force comparison? The new native input must equal job16602's exact template byte-for-byte after replacing its single `PAtom` token with `PModel`. No other scientific input change is permitted. The public calculator leaves the guess unspecified: this is a **PModel derived-recipe variant on ORCA6.1.1**, not exact historical ORCA6.0.0 reproduction.

Preserve the identical 256 physical IDs/coordinates, charge−2, UKS singlet, 1030 electrons, 515/515 occupations, La ECP46 and beta rotation `{514,515,20,1,1}`. Keep complete external orbital basis `def2-tzvpd.bas` SHA7389c56f0cff9e76098b0757a34d1abdd79c48e9249597cea074a0106a3462a6, built-in def2/J and external-file La ECP46 assignment, wB97M-V EnGrad, RIJCOSX, NoUseSym, DIIS/NOSOSCF, NormalConv plus Tight/300, THRESH1e−12, TCUT1e−13, DEFGRID3, population/output controls and NoAutostart. No prior GBW is read.

The unchanged EnGrad request caused native6.1.1 to select SCNL/libXC and disable AutoTRAH. Record and review that actual treatment; do not call it post-SCF VV10. No explicit SCNL keyword or grid override is added. PModel changes the initial orbitals and may change the converged SCF solution despite an unchanged final functional, basis and molecular state declaration.

## Finite stages and accounting

1. `prepare`: one CPU under Slurm; zero engine calls. Admit exact failed16602 result/config/source/input, inherited central geometry/state and source pins, external basis and exact one-token deck delta. Produce a fresh manifest and input evidence.
2. `run`: only after root binds the actual preparation and a decision with status `approved_single_pmodel_recipe_central`. Exactly one requested whole-molecule EnGrad endpoint through at most one native runner invocation; 24 actual MPI ranks, one worker, standard_expand default allocation memory/time and allocation-derived MaxCore. Fresh directories; no automatic retry, fallback or restart.

**Native PModel atomic-density subcalculations are within scope.** Prior15545 explicitly generated a La/ECP atomic density internally. Preserve all surviving task-local `*_atom*` files by hash in ATOMIC_GUESS_ARTIFACTS.json, along with native output and normal runner artifacts. One endpoint/runner invocation does not mean one total internal solver invocation; internal guess cost is included in native wall time and scheduler allocation totals. No separately requested atomic calculation is introduced.

No new MACE, ALPB, GFN, MM, DFT displacement, other-metal endpoint or optimization. No extra standalone preflight beyond the declared zero-call preparation.

## Admission and interpretation

Reuse unchanged native runtime/provenance, UKS/count/ECP46, normal-termination/SCF convergence, finite energy and full256 gradient gates. Final/native gradient energies must agree within1e−8Eh; gradient coordinates retain the existing1e−6Å comparison. Permit and record MODEL POTENTIAL; reject the old polarized-atomic-density guess on completed admission. Capture SCNL/libXC, rotation, basis/auxiliary/ECP, grid, density/gap/S² and residual evidence. Successful core admission remains `native_core_admitted_pending_method_review`; unknown details do not authorize rerunning chemistry.

Compare the saved full vacuum MACE and PBE0 forces with the new force and the previously frozen full256 direction. Preserve Cartesian-component RMS, maximum atom-vector norm, full vectors and projection units exactly as in the previous adapter. No cross-method absolute-energy work, physical-state/stability/spin-purity qualification, independently qualified force accuracy, shared-solvent validation, equilibrium or metal-preference claim. No automatic SP pair or further molecular release follows.

## Execution interface

`TARGET.py --stage prepare|run --config CONFIG --output FRESH_DIRECTORY`

CONFIG_PREPARE binds the implementation, this specification, original native/source assets, failed16602 lineage and pinned lanm_qmmm Python. `preparation_result` and `execution_decision` are null until actual evidence exists. Root owns launchers, actual run configuration, submission, monitoring and cost records. All scientific imports and calculations occur only inside Slurm; preparation of this source used only filesystem metadata/text and AST inspection.
