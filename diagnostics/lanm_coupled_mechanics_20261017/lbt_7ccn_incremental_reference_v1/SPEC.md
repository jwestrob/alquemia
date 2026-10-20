# One fixed incremental 7CCN-La reference

## Question and selection

Does the unchanged whole-peptide PBE0 plus cached fixed-FF ALPB model support the incremental lowering from the exact step80 endpoint to the exact terminal of the declared 12-step continuation? Root SCIENCE_DECISION.json (9503a623) freezes this question before execution. Selection uses the terminal of the complete declared scout, with no score-based intermediate selection.

One new native endpoint: POINT eb96df29, selected by actual continuation15254 RESULT c445cf4b. Its 12 accepted steps are physically admitted, still nonconverged. Reuse actual14880 after-La PBE0 endpoint RESULT816e2671 at old POINTcf59ef7a, already admitted in complete collector15173 fff95c2f. The model's saved incremental work is −4.226977666kcal/mol (MACE−1.705988079, ALPB−2.520989587). These values declare the comparison, not the reference outcome.

## Fixed model and state

Whole256 physical atoms, same exact IDs/topology/proton inventory, 7CCN Lu-conditioned source, La target index255, charge−2, physical/effective singlet. No caps or external field. Same PBE0-D4/def2-TZVP/def2-J/RIJCOSX/DefGrid3/VeryTightSCF/EnGrad/NoAutostart/PModel recipe and archived La basis/ECP46. Native input text must exactly equal the old14880 input template; all native assets, parser, runner, grid, SCF and memory-admission checks remain unchanged.

Reuse pinned c0cca6dd REFERENCE.inputs with the untouched original15173 config and REFERENCE.saved_run for actual14880. The inherited saved input checks admit four original cells without new molecular calls. Admit the new POINT's recipe against the old La recipe exactly, source/atom metadata, coordinates/XYZ/cache identity, MACE energy/full forces and the native ALPB receipt. Original whole-MACE gradient construction remains `-F_MACE*EV + g_ALPB*(HA/BOHR)`; native reference solvent gradients remain `(g_ALPB*HA)/BOHR` as produced by the original helper. Do not bitwise equate these differently associated representations.

The native admission function is copied unchanged from c0cca6dd, with the new one-cell manifest validator and explicit new-scope caveats. Reuse original whole-state, asset/deck, allocation, energy/analytic-gradient parser, canonical JSON and full-force-report helpers. No new energy engine or gradient algorithm.

## Stages, calls and resources

- `--stage preflight`: fresh output, Slurm one CPU, zero molecular calls/model loads. Re-admit old native endpoint, new saved POINT/components, both cached solvents, declared incremental model work and exact new single-task manifest. Future run bindings must be null.
- `--stage run`: exactly one fixed `new_La` PBE0 energy/analytic-gradient request. Requires the actual passed preflight and root's separately pinned execution decision with identical scientific binding. Twelve actual MPI ranks, one native worker, standard partition, standard_expand default wall time; allocation-derived memory, no memory/time override. No repeat/retry, optimization or other-metal task. Original native contract requires12 ranks; no24-rank branch is introduced.
- `--stage collect`: saved-only Slurm one CPU. Re-admit exact old/new native endpoints, preserve unavailable/failed new endpoint, and emit the incremental comparison only when both are available. Zero molecular calls.

No new MACE, ALPB, GFN, MM, geometry trials or optimizer steps in any stage. Inputs and output directories are immutable/fresh; errors preserve partial outputs, native receipt and both comparison slots. Root owns launchers, actual configs, submission, health monitoring and cost receipts. No stage launches another.

Execution decision schema: protocol `lbt_7ccn_incremental_reference_v1`, status `approved_single_fixed_reference`, exact `scientific_binding`, actual `preflight_result`, `permitted_execution_tasks:["new_La"]`, `new_point`, `old_reference`, and `mpi_ranks:12`. Config preflight_result/execution_decision/endpoint_result remain null until actual records exist.

## Report and interpretation

Take differences in native units before conversion: new-minus-old PBE0, MACE and cached ALPB; then composite works and MACE-minus-PBE0 work error. Report the pre-existing1kcal/mol material resolution and whether the reference composite lowers. This is neither a native admission threshold nor an uncertainty estimate. No absolute cross-method total comparison, affinity/preference or minimum claim.

Preserve all256 force vectors and model-minus-reference maximum/RMS atom-vector residuals at both endpoints. Report projections on the normalized exact incremental chord without alignment, geometry modification or alternate direction selection. These are descriptive: native DFT analytic gradients have not received independent numerical/accuracy qualification. Normal SCF termination does not override printed residual flags or establish stable state/electron binding. Retain vacuum-dianion, frozen-f/ECP, post-added unpolarized fixed-charge solvent and conditional terminal-chemistry limitations. A discrepancy does not identify a specific missing interaction. This single test authorizes no further relaxation or broader reference panel.
