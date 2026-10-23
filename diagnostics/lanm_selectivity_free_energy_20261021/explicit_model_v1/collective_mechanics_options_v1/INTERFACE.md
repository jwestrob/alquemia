# Staged trimer input contract

Preparation is complete and unsubmitted. `PREPARE.py` exposes no native execution method. Existing sources and admitted results are unchanged.

## Reproduce preparation in a fresh output directory

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/collective_mechanics_options_v1/PREPARE.py \
  --config diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/collective_mechanics_options_v1/CONFIG_PREPARE.json \
  --out workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/collective_mechanics_options_v1/preparation_v1
```

The displayed output already exists and is immutable; use its receipt rather than repeating it. The only executed work was source/metadata preparation and saved arithmetic, with zero native calls.

## Manifest fields for root's thin execution integration

`preparation_v1/MANIFEST.json` contains:

- `reference_cells`: ordered list of exactly16 cells, using the existing common-basis `REFERENCE.parse(text, cell, returncode)` input schema. Each has `template_input`, full `symbols/positions_A`, `ghost_mask`, `real_parent_atom_indices`, state/electron/AO/Aux-J/ECP inventories and parent source geometry. `common.deck(parent_panel, cell, 4, actual_maxcore)` renders the same recipe at the actual memory grant. Roles are `metal_two_water` and `two_water`; no original `definition()`/68-cell release gate can be reused unchanged.
- `reference_assets`, `orca`, `sources.common_panel/common_adapter`: exact reference assets, executable and helper pins. Copy or verify assets unchanged; do not replace the established run/release gate with the 68-cell gate.
- `model_rows`: exactly16 records keyed by new model ID. `input` pins XYZ/key/PRM; `source_input` retains the selected parent provenance. Schema matches qualified `RUN_v2.parse_moments(output, induced, row, induction)`, including `native_names`, `COM_A`, charge, site count and `energy_replay_kcal_mol=null` for these genuinely new fragments. Invoke pinned `model_analyze` with `system.xyz EM`; use native-announced induced filename and the same induction reader. No finite differences.
- `reused_reference`: exact admitted common-basis record plus the admitted common-origin quantum moment row for every old term used. Their original result/output pins are retained. `sources.common_collection/common_review/quantum_moments/quantum_moments_review` bind authority.
- `reused_model`: exact previous definition and actual native record for each old monomer/dimer/ion term. `sources.model_result/model_review` bind authority. No re-execution.
- `cycles`: eight records keyed by parent plus `near`/`far`. Each gives explicit `reference_terms` and `model_terms` coefficient dictionaries, geometry-only selection, atom/pair cancellation proof, and prospective model print bounds. Use these exact complete dictionaries; missing values invalidate the cycle.
- `selections`: original near/far O–O distance and M–O–O angle. A symmetry-suppressed far-pair dipole remains an observation, not a reason to replace the selected pair.

Root's native adapter must bind this manifest and a new narrow release for16 reference programs/four MPI ranks and at most16 model programs/four single-thread workers. The old78/68 adapters remain immutable; their parsers/renderers are reusable, their finite-inventory release gates are not. No launch configuration or approval is implied by `native_release=null`.

## Interpretation and next implementation

The report recommends a metal-pair-scoped distinction between direct driving fields and mutual induced response, implemented through one consistent energy functional and analytic force path. It is a candidate architecture, not a fitted or built model. This trimer panel measures where its prospective constraints arise. The test can run independently of the pending correlated-reference comparison because it diagnoses the provisional PBE0 surface; it cannot promote that surface to an accurate reference.
