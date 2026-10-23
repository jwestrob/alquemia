# Ready: matched DFT counterpoise comparison

Five fixed cells, four MPI ranks each; root owns submissions. Native pilot and
parser evidence are in SPEC.md and validation_v1/VALIDATION.json. No native
calls were made during preparation. All 22 parser controls passed, including
both real complex outputs, all five roles and 15 negative controls.

- Panel: `workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dft_water_cp_v1/prepared_v3/PANEL.json`
- Panel SHA256: `3efe421275c6e9833eea757b6b0990e04195283a571326f65c40708ece1f1d0e`
- REFERENCE.py SHA256: `b1f35c46cf4512b6088d5a81183cc6ab1532fd997fee05ab5fd904b3c4595089`
- Validation SHA256: `705eaa22594584ea08461112ecf51b99b87f2ab132da7b5c43319f87116b1157`
- Cells: `La__ion_ghost_water`, `La__water_ghost_ion`, `Dy__ion_ghost_water`, `Dy__water_ghost_ion`, `real_water_D4`.

Pass these arguments after the owned RUN.sbatch wrapper:

    /absolute/path/REFERENCE.py run --panel /absolute/path/PANEL.json --cell CELL --out /new/absolute/execution

Use unique outputs. Native checks must establish ghost basis/charge/ECP and
no-D4 behavior for four ghost roles, and three real atoms/43 basis functions
plus actual D4 for water. Collection uses the complete canonical collection:

    REFERENCE.py collect --panel PANEL --canonical CANONICAL_COLLECTION --out NEW_RESULT RESULT1 RESULT2 RESULT3 RESULT4 RESULT5

Original complex calculations are reused, with printed electronic and D4
energies separated and closure checked. The isolated-ion D4 term is an explicit
empty-pair/triple sum, not a native result. The real-water D4 term is retained
in each individual interaction even though it cancels between metals.

Preparation attempts v1/v2 failed before molecular calls on legacy result-schema
keys. Their files and source snapshots remain preserved. Prepared v3 succeeds.
The reference/scientific settings, source recipes and previous calculations
were unchanged by those preparation repairs.
