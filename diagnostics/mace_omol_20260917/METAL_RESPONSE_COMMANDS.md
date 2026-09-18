# Full metal response operations

Production baseline and previous two-coordinate experiment remain unchanged.
Read METAL_RESPONSE_PLAN.md for fixed inputs, inventories and acceptance checks.
Current preparation is workspaces/mace_metal_response_20260918/prepared_v2.
V1 was not executed; V2 freezes its preparation source for later development.

From the repository root, inspect/collect an existing job without new inference:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_metal_response.py validate --manifest workspaces/mace_metal_response_20260918/prepared_v2/core_GGR_extended/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_metal_response.py collect --manifest workspaces/mace_metal_response_20260918/prepared_v2/core_GGR_extended/manifest.json --output workspaces/mace_metal_response_20260918/prepared_v2/core_GGR_extended/recollection_01.json
```

All seven actual submission.json files contain exact commands, manifests and
job IDs. Existing run_pilot.sbatch invokes the frozen implementation snapshot;
one A5000,16CPU,64474MiB per job. Reuse completed accepted attempts on recovery;
never launch a second executor for a live manifest or erase a failed attempt.
The original center energy/gradient/short/GB records remain reused and pinned.

After core collection, prepare the matched solvent tasks (fresh output only):

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_metal_response.py prepare-gb --collection workspaces/mace_metal_response_20260918/prepared_v2/core_GGR_extended/collection_job_1201352.json --solver-validation workspaces/mace_gb_20260916/pilot_v2/collection_job_1200700.json --output workspaces/mace_metal_response_20260918/gb_GGR_extended_v1
```

The other three core manifests use jobs1201353,1201354,1201355. Whole short
manifests use1201356,1201357,1201358. Every job has72 tasks; no DFT hidden in them.
Use the pinned solvent environment workspaces/mace_gb_20260916/software_v2/venv/bin/python
for GB execution. The worker records native OBC2 reaction energy, never adds
its fixed-charge Cartesian force to a learned-charge total-gradient claim.
Curvature comes from newly evaluated charge distributions at every geometry.

Assessment takes explicit preparation, JSON index with a `collections` list
of completed core/short/GB collection paths, and a fresh output directory.
It rejects missing/failed tasks. Minima are predictions until a separate native
analytic DFT+short validation passes. `mace_metal_minimum.py prepare` takes the
assessment result and fresh output, preserving all excluded endpoints.
No PQQ threshold, entropy, or automatic score promotion exists here.

## Completed grid and independent DFT checks

All792 declared cheap evaluations are complete (288core,216whole-short,288GB).
The explicit eleven-collection index is collections_v1.json. Read-only replay:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_metal_response.py assess --prepared workspaces/mace_metal_response_20260918/prepared_v2/preparation.json --collections workspaces/mace_metal_response_20260918/collections_v1.json --output workspaces/mace_metal_response_20260918/assessment_replay_v1
```

Actual assessment_v1 finds all eight combined matrices positive and converged.
All four alpha optima exceed0.20A and remain unavailable. All four GGR optima
are eligible; minimum_v2 contains exactly4native analyticDFT +8short tasks.
minimum_v1 was not executed: GPU-environment preflight found CPU/GPU NumPy
roundoff up to1.14e-13 in recomputed eigensystem metadata. V2 allows1e-10
floating comparison while retaining exact statuses and original predicted
coordinates. Physical acceptance tolerances are unchanged. Three real-grid/
input/malformed-result tests pass0.852s; frozen GPU preflight passes.

Actual quantum job1201370 and short job1201371 use existing runners. Their
submission.json files contain complete commands and hashes. After completion:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_metal_minimum.py report --prepared workspaces/mace_metal_response_20260918/minimum_v2/preparation.json --dft workspaces/mace_metal_response_20260918/minimum_v2/quantum/collection_job_1201370.json --short workspaces/mace_metal_response_20260918/minimum_v2/short/collection_job_1201371.json --output workspaces/mace_metal_response_20260918/minimum_report_v1
```

Endpoint-validated and fully qualified responses are separate. The latter
remain null unless the required GGR partition check also passes. Unvalidated
actual changes remain available under explicitly exploratory fields.

Both validation jobs are now complete. Primary immutable result: minimum_report_v2/result.json; V2 pins the reporter source locally and preserves identical scientific values. Older partial/final reports and their original reporter copies remain available.
