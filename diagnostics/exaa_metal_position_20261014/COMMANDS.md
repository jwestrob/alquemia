# Exact operations for the current finite test

Run from the repository root. Preparation4129 is complete and scoring4130 is
already submitted; these are reproducibility/recovery commands, not a request
to duplicate either job. Inspect CURRENT.md, EXECUTION.json and actual scheduler
state first. Immutable existing output directories deliberately reject reruns.

The exact preparation submission, manifest hash and wake command are in:
`workspaces/exaa_metal_position_20261014/preparation_v1/SUBMISSION.json`.
The exact scorer submission, config hash and wake command are in:
`workspaces/exaa_metal_position_20261014/scoring_v1/SUBMISSION.json`.

Read actual terminal products without rerunning molecular work:

```bash
cat workspaces/exaa_metal_position_20261014/preparation_v1/execution/PREPARATION.json
cat workspaces/exaa_metal_position_20261014/scoring_v1/execution/EXECUTION.json
cat workspaces/exaa_metal_position_20261014/scoring_v1/execution/FINAL_COLLECTION.json
```

All scientific preparation, checking and scoring belongs on Slurm. The preparation
manifest pins the exact source inventory and frozen physical helpers. The scorer
config pins preparation, model, runtime, archived qualification and original
energies. No activation, source editing, fresh download or hidden defaults are
required by the actual recorded submissions. Missing/partial outputs are retained
as such and are not a successful fallback to old energies.

## Current solvent continuation4132

The completed MACE and exact original controls are bound in
`workspaces/exaa_metal_position_20261014/solvent_v1/CONFIG.json`.
The actual submitted argv and completion-wake command/PID are retained in its
`SUBMISSION.json`. Do not repeat metadata prepare or execute into this workspace.
Actual internal preflight passed before new SCC work.

Read the terminal result without rerunning calculations:

```bash
cat /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/exaa_metal_position_20261014/solvent_v1/execution/FINAL_COLLECTION.json
```

The collector reports original/transferred per-metal donor works, their MACE and
solvent components, D and delta_D, and separate complete two/four-geometry row
minima. No class bands or calibration are applied. Missing required cells stay
unavailable. The unrounded records and hashes are authoritative.
