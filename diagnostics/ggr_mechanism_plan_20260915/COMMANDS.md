# Reproduce and inspect the approved GGR investigation

Repository root:
`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs`.
The script below records that root and the existing Python executable explicitly.
It recollects pinned outputs, checks receipts, produces comparisons/components,
records actual Slurm costs, exports a derived benchmark ledger and renders four
PNG/PDF figures. **It launches no quantum calculations.** Its output directory
must be new; original results are never overwritten.

```bash
bash /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/ggr_mechanism_plan_20260915/rebuild_report.sh \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ggr_mechanism_20260915/report_recheck_v1
```

After that path exists, pass another new absolute workspace directory. Collection
preserves missing/failed statuses; it cannot turn an unfinished job into a result.

## Frozen tasks and actual submissions

| Stage | Manifest under `workspaces/ggr_mechanism_20260915/` | Job | Collector |
|---|---|---|---|
| A | `stage_a_tasks_v1/manifest.json` | 1199770 | `ggr_workflow.py collect` |
| B | `stage_b_tasks_v1/manifest.json` | 1199802 | `ggr_workflow.py collect` |
| C | `stage_c_tasks_v1/manifest.json` | 1199805 | `ggr_sensitivity.py collect` |

The exact executed submission argument arrays and hashes are in
`SUBMISSION_A_1199770.json`, `SUBMISSION_B_1199802.json`, and
`SUBMISSION_C_1199805.json`. Each task directory retains its prepared input,
coordinates, output, runtime input and hash-linked execution receipt. Gradient
centers also retain the actual `.engrad`. Implementation snapshots accompany
the manifests; collections additionally record the actual collector/adapter.

No resubmission is needed to rebuild a report. For an interrupted agreed task,
the existing execution runner accepts the same manifest and checks completed
receipts before reusing them. It refuses to overwrite incomplete/unverified
attempts; preserve such an attempt before any unchanged-science technical
recovery. Do not reuse this task list to add a new representation or structure.

Preparation commands and exact frozen inventories are in
[STAGE_A_PREPARATION.md](STAGE_A_PREPARATION.md) and
[STAGE_B_PREPARATION.md](STAGE_B_PREPARATION.md). Gradient collection and the
conditional half-step rules are in [GRADIENT_INTEGRATION.md](GRADIENT_INTEGRATION.md)
and the [approved plan](PLAN.md). The half-step command selects blocks from
verified failed numerical checks, never from a score's sign or missing output.

See [TESTS.md](TESTS.md) for executed software/integration checks, and
[REPORT.md](REPORT.md) for scientific interpretation and the next proposed
scientific decision. The latter is not an additional execution authorization.
