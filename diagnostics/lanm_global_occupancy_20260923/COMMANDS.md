# Exact pilot operations

Run from the repository root. Jobs1213018 and1213040 have already been submitted;
do not repeat their submission commands. Their actual finite preparation and
implementation are pinned in `workspaces/lanm_global_occupancy_20260923/`.

Current queue and local checks:

```bash
squeue -j 1213018,1213040 -o '%.18i %.25j %.8T %.10M %R'
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/lanm_global_occupancy_20260923/test_scoring.py
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/lanm_global_occupancy_20260923/scoring_v2/implementation/affordable_workflow.py dry-run --manifest workspaces/lanm_global_occupancy_20260923/native_feasibility_v2/manifest.json
```

After the first job finishes, this read-only operation verifies actual feasibility
without launching anything:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/lanm_global_occupancy_20260923/continue_pilot.py gate
```

The dependent job executes the remaining eight prepared systems automatically if
that gate passes. The final native CPU stage uses the existing manifest executor;
its exact job ID is recorded in `workspaces/lanm_global_occupancy_20260923/NATIVE_SUBMISSION.json`.
It writes the final JSON/report/vault/email automatically, including failed cells.

For an independent collection after those artifacts exist, with no email and
no new calculations, choose the distinct output below (refuses overwrite):

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/lanm_global_occupancy_20260923/scoring_v2/implementation/lanm_global_occupancy.py collect --manifest workspaces/lanm_global_occupancy_20260923/scoring_v2/manifest.json --mace-result workspaces/lanm_global_occupancy_20260923/mace_merged_v2/result.json --native-manifest workspaces/lanm_global_occupancy_20260923/native_matrix_v2/manifest.json --output workspaces/lanm_global_occupancy_20260923/independent_collection_v2.json
```

Use `--help` on the pinned `lanm_global_occupancy.py` for explicit preparation,
MACE and native-preparation arguments. Do not reuse an output directory or run
fresh chemistry merely to collect an existing result. Actual integration remains
pending until scheduler execution produces molecular receipts.
