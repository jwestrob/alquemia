# Bounded metal-response operations

Research protocol `DFT_anchored_MACE_GB_bounded_metal_response_v1`.
Read the frozen BOUNDED_METAL_RESPONSE_PLAN.md. Production remains unchanged.
Commands run from the repository root. Output directories must be fresh.

## Prepare and inspect (no scientific execution)

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_bounded_response.py prepare --assessment workspaces/mace_metal_response_20260918/assessment_v1/result.json --reuse-result workspaces/mace_metal_response_20260918/minimum_report_v2/result.json --agreement diagnostics/mace_omol_20260917/BOUNDED_METAL_RESPONSE_PLAN.md --output workspaces/mace_bounded_response_20260918/prepared_replay_v1
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_bounded_response_20260918/prepared_v3/short/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_bounded_response_20260918/prepared_v3/short/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/affordable_workflow.py dry-run --manifest workspaces/mace_bounded_response_20260918/prepared_v3/quantum/manifest.json
```

Actual preparation is `prepared_v3`: four alpha native analytic DFT tasks and
eight short MACE tasks. All four GGR states reuse exact archived displacements,
energies and gradients. No new curvature grids. All donor memberships pass.
V1/V2 were unexecuted: GPU preflight exposed a preparation-only import absent
from the inference environment. CPU preparation/report fully validate archived
DFT receipts; GPU preflight verifies their pinned artifacts without importing
preparation-only dependencies. V3 changes no scientific inputs or predictions.

## Execute and collect

Existing runners submitted jobs1201383 (quantum) and1201384 (short). Each
`prepared_v3/{quantum,short}/submission.json` contains the exact executable argv,
manifest hash and scheduler response. Quantum uses four16-rank workers on64CPU;
short uses one A5000,16CPU,64474MiB. Do not start a duplicate executor for a live
manifest. Successful accepted attempts remain reusable; failed attempts persist.

Read-only recollection after execution:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_mechanics.py collect-dft --manifest workspaces/mace_bounded_response_20260918/prepared_v3/quantum/manifest.json --output workspaces/mace_bounded_response_20260918/prepared_v3/quantum/recollection_v1.json
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_bounded_response_20260918/prepared_v3/short/implementation/mace_hybrid.py collect --manifest workspaces/mace_bounded_response_20260918/prepared_v3/short/manifest.json --output workspaces/mace_bounded_response_20260918/prepared_v3/short/recollection_v1.json
```

## Compare and report

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_bounded_response.py report --prepared workspaces/mace_bounded_response_20260918/prepared_v3/preparation.json --dft workspaces/mace_bounded_response_20260918/prepared_v3/quantum/collection_job_1201383.json --short workspaces/mace_bounded_response_20260918/prepared_v3/short/collection_job_1201384.json --output workspaces/mace_bounded_response_20260918/report_replay_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_bounded_response.py -v
```

Reports retain actual exploratory and qualified fields separately. Qualification
requires energy/gradient checks and the GGR partition check. Boundary tangent
and radial conditions describe the fixed sphere, not unconstrained equilibrium.
Both alpha-minus-GGR margins and directions are reported without calibration.
No inherited PQQ bands, entropy correction, production score or default change.
