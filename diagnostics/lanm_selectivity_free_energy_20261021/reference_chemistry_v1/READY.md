# Eight native pilot cells ready

Root reviewed the finite panel and runner. No molecular calculation was run by
the preparation agent. The root owns all Slurm launches and terminal admission.

Use `workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/prepared_v2/PANEL.json`.
Its SHA256 is `d5929526cf70ac957c847a5936ed229e9edf1322764b0cc287a59292e01204f9`.
The first `prepared/` panel is retained as the early physical-coordinate handoff;
its adapter pin predates the native NEL/Mult parser correction. Coordinates are
unchanged. Use only prepared_v2 with the final native runner.

The eight IDs are in `pilot_cell_ids`. For each singleton, call this directory's
`RUN.sbatch` with arguments:

```
/absolute/path/to/REFERENCE.py run --panel /absolute/path/to/prepared_v2/PANEL.json --cell CELL_ID --out /absolute/path/to/new_cell_execution_directory
```

The script requests four MPI ranks, standard_expand and scheduler-default memory
and time. The runner derives MaxCore from the actual grant and requires Slurm.
Root must arm actual completion notifications. Native execution creates START,
input, output, error, gradient and immutable RESULT records per cell.

Local preparation checked28 exact paired physical fixtures and rendered56 native
decks with zero molecular calls. Native parser replay on an actual archived Dy
output passes. Five negative controls reject wrong electron count, multiplicity,
core treatment, coordinates and failed process status. See
`PREPARATION_VALIDATION_V2.json` and `NATIVE_PARSER_CONTROLS.json`.

After admitted pilot results, `collect` accepts their RESULT paths and reports
the balanced four-cell acetate-versus-water contrast and the difference between
electronic recipes. Missing cells remain missing; this is a partial56-cell panel.
Full per-atom forces remain in each result. No automatic native retry, model
promotion or claim of reference accuracy follows from process completion.
