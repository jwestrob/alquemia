# Stage A preparation completed

Five source-pinned pairs are ready for the parent workflow (10 endpoints).
No ORCA calculation was submitted by the preparation module.

| Preparation | Atoms | Source displacement (Å) | Ligand charge |
|---|---:|---:|---:|
| ggr_1glg_nma | 58 | 0.0 | -3 |
| aequorin_1sl8_EF3_nma | 49 | 0.0 | -3 |
| alacta_1f6s_nma | 52 | 0.0 | -3 |
| alacta_6ip9_nma | 55 | 0.0 | -3 |
| ggr_1glg_connected | 111 | 0.0 | -3 |

All parent v3 source/atom/cap coordinates were reconstructed exactly at their serialized precision before extending each representation. Paired coordinates, source hydrogen placement, explicit waters, formal charges and donor records remain fixed. Graph overlaps merge; GGR’s connected segment contains seven caps and no Lys137 side-chain atoms. Endpoint charges remain La 0/Ca −1, singlets.

The two experimental protocol IDs are `generic_peptide_alpha_caps_native_r2scan3c_dev_v1` and `ggr_connected_segment_native_r2scan3c_dev_v1`. Neither inherits a decision cutoff or compatible absolute reference. Default production code was not changed.

Eight real-fixture tests passed: all pair invariants, actual NMA anchor/cap topology, connected overlap/charge bookkeeping, deterministic preparation/cache identity, snapshot integrity, failure on a corrupted real missing-anchor input, rejection of the charge-changing aequorin connected segment, and unchanged input/XYZ hashes for all 50 historical PQQ endpoints. These tests ran no electronic-structure executable.

The only observed warning was an existing `ElementTree` truth-value deprecation in the borrowed `affordable_peptide.py`; that shared file was not edited.

## Runnable operations

Replay to a new output directory (existing outputs are immutable):

```bash
python scripts/ggr_preparation.py batch \
  --config diagnostics/ggr_mechanism_plan_20260915/STAGE_A_PREPARATION_CONFIG.json \
  --output workspaces/ggr_mechanism_20260915/stage_a_prepared_replay
```

For one prepared v3 parent, use `prepare --repair-manifest PATH --policy alpha_caps --output NEWDIR`. The same API is `prepare(repair_manifest, output, policy="alpha_caps")`; `connected_segment` is restricted to the approved topology/sequence-defined GGR chemistry.

Each directory contains `preparation_manifest.json`, both XYZ/input files and immutable implementation copies. `atom_graph.source_to_qm` retains the existing convention: ligand index 1 is full XYZ index 1 after prepending metal index 0.
