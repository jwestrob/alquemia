# One reference at the model-converged 7CCN–La endpoint

## Question and fixed selection

Does unchanged whole-256-atom PBE0 plus the same cached fixed-FF ALPB support the **step 92 to step 725** model energy lowering? What reference residual forces remain at the model-force-converged terminal?

Root SCIENCE_DECISION `27377db9c2064650094e9f6106a106b57509a3701fa233cc7dfda5065e1be83b` declares exactly one new endpoint. It is the original-force-criteria terminal selected by actual continuation 15525, RESULT `95d870eb`, POINT `6605f43c`, checkpoint `0461e30d`. No intermediate or outcome-favorable geometry is selected. The stored model maximum/RMS forces are 0.0750547504/0.0296109468 kcal/mol/Å; they establish model convergence only.

Reuse actual 15321 PBE0 energy/gradient at step 92, RESULT `3a173d81`, admitted by collector 15456 `5908cb2d`. Its POINT is `eb96df29`. Reuse both saved MACE energies/full forces and both fixed-FF ALPB receipts (step92 `37e90528`, step725 `9a3b7c6e`).

## Unchanged chemistry and native protocol

Same 256 physical IDs, topology/protons, Lu-conditioned 7CCN source, La target index255, charge−2, physical/effective singlet, no caps or external field. Preserve the exact previous PBE0-D4/def2-TZVP/def2-J/RIJCOSX/DefGrid3/VeryTightSCF/EnGrad/NoAutostart/PModel input and archived La basis/ECP46. The generated input text must equal the actual 15321 native template byte for byte, including all original opaque settings.

Use unchanged incremental adapter `bf75f502` with its exact actual collection config `2fef7587` for saved input and 15321 native admission. Its former `new_La` cell becomes this comparison's `old_La` without modifying the original native row, task ID, receipts or artifact pins. Retain its original `c0cca6dd` helper for task signatures, runtime and force reports. The new adapter's native `admit` body is copied unchanged; parser, asset checks, runner, canonical row handling, SCF semantics and allocation-derived memory checks remain fixed.

Admit the new POINT's exact source/atom metadata, recipe, state, coordinates/XYZ/cache identity, raw MACE state/energy/full forces and ALPB provenance. Preserve original composite-gradient arithmetic `-F_MACE*EV + g_ALPB*(HA/BOHR)` separately from the native reference solvent representation `(g_ALPB*HA)/BOHR`. No tolerance or representation substitution is introduced.

## Stages, calls and gates

- **Preflight:** one Slurm CPU, fresh output. Re-admit the old reference and new saved components, both cached solvents, exact terminal/checkpoint/source lineage, and the one-task native manifest. Zero molecular calls or model loads. Future execution bindings remain null.
- **Run:** exactly one fixed `new_La` PBE0 energy/analytic-gradient call, after actual matching preflight and root's pinned execution decision. Twelve actual MPI ranks, one native worker, standard partition and `standard_expand` defaults. Allocation-derived MaxCore; no memory/time overrides or resource-contract changes. No retries.
- **Collect:** saved-only one-CPU Slurm job after the actual endpoint exists. Re-admit both native endpoints and produce the comparison only if both are available. Preserve failed/unavailable slots and native evidence. Zero molecular calls.

No MACE, ALPB, GFN or MM calls; no optimization, geometry trials, other metal/source or other reference cells. All numerical imports/replay run under Slurm. Source/JSON/hash/AST inspection alone is allowed on login. Root owns launchers, actual run/collection configs, submissions, watchers and scheduler cost records. No stage launches another.

Execution decision schema remains: this protocol; status `approved_single_fixed_reference`; exact `scientific_binding`, actual `preflight_result`, `permitted_execution_tasks:["new_La"]`, `new_point`, `old_reference`, and `mpi_ranks:12`. `preflight_result`, `execution_decision` and `endpoint_result` are null in CONFIG_PREFLIGHT.

## Comparison and limits

Take **native endpoint differences before unit conversion** for PBE0, MACE and ALPB. Job15525's `work_from_checkpoint` spans only step172→725 and is not the step92→725 comparison. Report both composite works, model-minus-reference work error, and the existing 1 kcal/mol reporting resolution. The resolution is neither a native admission gate nor an uncertainty estimate.

Preserve all256 force vectors, endpoint maximum/RMS forces, model-minus-reference force-error metrics and projections on the normalized step92→725 full-coordinate chord. No alignment or direction selection is introduced. Forces remain descriptive and independently unqualified; model convergence is not reference stationarity or a physical minimum.

Shared cached ALPB cannot independently validate solvent physics. Retain conditional termini/protonation, vacuum-dianion electron-binding uncertainty, post-added fixed-charge solvent without density polarization, frozen-f/ECP assumptions and all printed SCF residual warnings. Native normal termination does not erase those flags. No affinity, metal preference, production promotion or automatic follow-up is authorized by this single comparison.
