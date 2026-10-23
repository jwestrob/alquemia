# Root integration

`PROBE.py prepare --config CONFIG_PREPARE.json --output FRESH_DIRECTORY` admits the exact failed-science terminal selection and its passed independent raw review without molecular calls. It writes RESULT.json, MANIFEST.json, compact DATA.json and 40 typed input directories. Previous preparation/source versions remain preserved.

The manifest's `data` pin resolves to DATA keys `rows`, `baseline_model`, `baseline_score`, `baseline_parameters`, `selected_normalized_parameters`, `family`, `six_column_Jacobian`, `six_column_actual_spans`, `derivative_keys`, `admitted_native_authorities`, `provenance`. Rows retain the original ten geometry/input/body/reference records. J6 is exactly reconstructed at the selected point from all twelve saved parameter derivative vectors.

Manifest `probes` order is La_minus5, La_minus1, La_plus1, La_plus5, then the corresponding Dy IDs. Each entry has `metal`, `alpha` (Å³), `eta` (actual relative alpha displacement), declared `fraction`, five `inputs` records and five `reused_cells` GIDs. Inputs contain exact XYZ/key/parameter pins, baseline input pins and one-token alpha-change evidence. No repeated large parent-result read occurs inside workers.

`PROBE.py run --config CONFIG_RUN.json --output FRESH_DIRECTORY` requires an eight-CPU allocation, all five thread variables equal to one, actual preparation and root release. The release fields are:

```
status: approved_fixed_metal_alpha_diagnostic
implementation: exact PROBE.py pin
preparation_result: exact preparation RESULT.json pin
scientific_config_sha256: fingerprint excluding root_release/preparation_result
maximum_new_native_programs: 80
maximum_new_physical_cells: 40
workers: 8
candidate_selection_authorized: false
```

The run config is the prepare config with `preparation_result` and `root_release` pins filled. Runtime RESULT contains `manifest` and `preparation_result` pins, `records[probe_id]` with status/metal/alpha/eta/fraction, all ten `cells` slots, complete `model` and unchanged `score` including `weighted_residuals`. Failed or absent slots remain explicit. Native files/receipts are in `PROBE_ID/GID/`; changed-cell models retain actual runtime-input/native-implementation provenance, and reused cells keep their original producers.

Terminal `metal_alpha_diagnostic_complete` means all eight native matrices and scores are available; it does not mean scientific controls passed. `failed_or_unavailable` preserves global stop and error. Cost fields are native invocation count, wall, parent CPU and aggregate child user/system CPU. Root owns separate two-scale column analysis, release, launch and any scientific follow-up.

The executable ROOT_DECLARATION/SPEC use independently computed relative rcond=1e−6 cutoffs for the original six and augmented eight matrices. This supersedes the fixed-base-cutoff suggestion in the earlier mechanics proposal; no additional algebra variant is requested. Report the implied unexecuted directions as diagnostics only.
