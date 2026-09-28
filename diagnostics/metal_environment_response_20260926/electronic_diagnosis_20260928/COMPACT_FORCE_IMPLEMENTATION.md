# Compact gradient diagnostic implementation

Prepared before inspecting compact gradient outcomes. Both existing workers1220308/1220312 remained running when this interface was completed; no collector was duplicated, no molecular call or submission occurred.

The frozen COMPACT_FORCE_DESIGN.json records each paired geometry hash, selected real oxygen/index/distance/direction and implementation hash. `scripts/metal_environment_compact_forces.py` reads an existing collection, verifies each underlying engrad hash, full atomic order, coordinates, energy and unit-converted gradient agreement. Missing/failed rows are unavailable; a partial pair never supplies a differential load. The saved Dy origin is explicitly reused.

Three real-fixture parser/algebra tests pass: source geometry/translation direction, actual Dy gradient sign/projection/translation sum, and unavailable or deliberately corrupted saved gradient record. These test analysis correctness, not new molecular force qualification. No new tolerances or displacement optimization are introduced.

Run after the existing final collection is available, from repository root:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/metal_environment_compact_forces.py analyze --design diagnostics/metal_environment_response_20260926/electronic_diagnosis_20260928/COMPACT_FORCE_DESIGN.json --collection workspaces/metal_environment_response_20260926/lady_compact_exchange_v1/FINAL_COLLECTION.json --output diagnostics/metal_environment_response_20260926/electronic_diagnosis_20260928/COMPACT_FORCE_RESULT.json
```

The JSON and adjacent Markdown report are exclusive-create outputs. Use a new explicit filename for a later collection rather than overwriting an analysis. Root owns scientific interpretation and promotion. Large individual loads suggest local strain on the specified conditional surface; neither loads nor their difference measure preference, relaxation energy, populations or full protein mechanics.
