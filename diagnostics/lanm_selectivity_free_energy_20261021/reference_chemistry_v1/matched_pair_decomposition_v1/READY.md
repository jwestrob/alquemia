# Literal 36-cell own-basis reference panel

Preparation and 54 targeted saved controls pass with zero native calls. `READY.json` pins the source, configuration, actual `prepared_v1/PANEL.json`, preflight, controls, wrapper and launch/release templates. Independent review and root release remain required before execution.

The panel contains 34 dimers, each an unchanged atom subset of the consumed W8/W9 parent, plus two real bare ions at the common metal origin. No symmetry replication, rotation, recentering or ghost basis is used. Every cell retains PBE0-D4/def2-TZVP/def2-J/RIJCOSX/DefGrid3, the archived real-metal orbital/auxiliary basis and ECP46/55, effective charge +3/singlet and original SCF settings. Dimers have 90 AO/18 electrons; own ions 47 AO/8 electrons. No EnGrad or force target is requested. Bare-ion D4 is evaluated by the declared real-atom recipe and its printed correction must be exactly zero; this is a future native admission check, not synthetic evidence that the atom calculation has run.

## Common geometry interface

`PANEL.cells` lists `role=dimer` (34) and `role=bare_ion` (2). Each contains `geometry`, its ordinary XYZ pin, `parent_geometry_id`, one-based `water_index`, zero-based `parent_atom_indices`, full native inventory and the frozen input-template pin. Each dimer maps `[0,3i−2,3i−1,3i]` from its original parent. The model author uses exactly these rows to construct the typed four-atom inputs.

Top-level `parent_geometries`, `saved_complexes` and `baseline_model_energies` retain the four consumed cluster authorities. The CP collection is context only; none of its ghost-fragment energies enters these own-basis cells. Original targets remain unchanged.

## Root launch and collection interface

`ROOT_RELEASE_TEMPLATE.json` has null status/review. Root binds the exact independent review and sets `approved_matched_pair_36_cells`, preserving the implementation, panel, preflight, full 36-cell inventory, four ranks and execution root. `LAUNCH_TEMPLATE.json` contains every exact `sbatch` argv and expected native output/result path. Create task parents before submitting; leave each `execution` directory absent so the adapter can create it exclusively. Existing attempts are rejected; no automatic retry occurs.

The unchanged four-MPI wrapper requests standard_expand with scheduler-default memory/time. MaxCore is 75% of actual memory divided by four; the template uses the established 16,000 MiB grant/MaxCore3000. Root schedules within current shared limits and arms completion watchers.

Public `load(panel)`, `deck(panel,cell,ranks,maxcore)` and `parse(output_text,cell,returncode)` support root's separate saved collector. Each RESULT records exact cell/panel/source/release/START/input/output/stderr, native admission fields, electronic/D4/final energy, SCF/stage metadata, return code and elapsed×ranks allocation cost. Failed attempts remain explicit results; unavailable entries cannot complete the 36-term decomposition.

The preflight records the executed preparation command and 0.952 s wall time; the single completed saved control set records its command and 1.902 s wall time, one local thread. Controls cover all 36 literal mappings/templates, four full per-parent denominators, both existing real-dimer native parser fixtures, the isolated D4 guard, malformed inventories and unreleased/duplicate attempt gates. New native own-ion and dimer admission remains pending. No collector or molecular execution was performed here.
