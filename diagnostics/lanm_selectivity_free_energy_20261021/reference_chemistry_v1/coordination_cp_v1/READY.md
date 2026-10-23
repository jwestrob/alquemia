# Eight-cell coordination counterpoise preparation

Eight electronic single-point decks are prepared and checked; **no molecular call has run**. Use the exact artifact pins in `READY.json`. `prepared_v4` is the final preparation; prior drafts and zero-call schema failures remain preserved.

## Scope and admission

For each consumed La/Dy W8/W9 complex, run the real ion with the whole shell ghosted and the intact real shell with the metal ghosted. The existing donor counterpoise renderer/parser supplies the same PBE0, basis, auxiliary basis, ECP and ghost conventions. The parser adds retained native final-integration metadata. D4 and EnGrad are absent from all eight fragment decks; the four original real-complex electronic/D4 energies remain separately pinned.

Prepared templates assume four MPI ranks and MaxCore3000. At execution, MaxCore is 75% of the actual Slurm memory grant divided by four; no memory or time override is requested. `RUN.sbatch` is byte-identical to the existing donor wrapper: standard partition, standard_expand, four one-thread MPI ranks, one node, no requeue.

Preparation replays all four original raw energy/gradient results, retains exact consumed geometry bytes despite the historical W8 path alias, and proves own-ion and own-whole-shell cancellation conditions. The targeted saved controls replay four real cluster outputs and four prior dimer ghost outputs, check all eight rendered coordinates, and reject wrong native/state/basis/ECP/D4/energy/termination and release/attempt evidence. These controls do not pre-admit the future cluster ghost outputs.

## Root launch interface

`ROOT_RELEASE_TEMPLATE.json` is deliberately unreleased. Root must bind its independent review pin and set status to `approved_coordination_cp_eight_cells`, retaining the exact implementation/panel/preparation, eight cell IDs, four ranks and execution root. Save the concrete release at the path given in `LAUNCH_TEMPLATE.json`.

`LAUNCH_TEMPLATE.json` supplies all eight exact `sbatch` argv arrays and expected native output/result paths. Create each task parent directory before submission; do not create its `execution` directory. The native adapter creates that directory exclusively, admits only the released cell path, and performs one attempt. An existing attempt is rejected. Root owns submissions and completion watchers.

Per-cell CLI:

```text
/usr/bin/python3 REFERENCE.py run --panel /absolute/prepared_v4/PANEL.json --release /absolute/run_v1/ROOT_RELEASE.json --cell CELL_ID --out /absolute/run_v1/tasks/CELL_ID/execution
```

Each result preserves actual input/output/stderr pins, native checks and stage metadata, release/source/panel/cell authority, return code, elapsed time and elapsed×four allocated CPU-seconds. Failure remains a result; there is no retry. The separate saved collector is owned by Arendt and uses `load(panel)` and `parse(raw_text, cell, returncode)`.

## Preserved preparation history

The first preparation stopped before native work because the W8 consumed reference dictionary lacks the W9-only optional `raw_result` key. Its authoritative `reference_result` remains mandatory and raw replayed. The second stopped on the historical W8 XYZ path alias; the repair verifies both files, identical hashes and every nonpath geometry field. A saved-control fixture-path typo pointed to a nonexistent old DFT preparation; its failed receipt/source is preserved. No scientific setting or acceptance criterion changed.

The final actual preparation command and timings are in `prepared_v4/PREFLIGHT.json`; the one completed control set and command are in `saved_controls_v2/RESULT.json`. Do not rerun preparation over existing paths. Original targets and the failed joint model remain unchanged.
