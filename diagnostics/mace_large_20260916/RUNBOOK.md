# Large analytic MACE pilot

Standing authorization: [AGREEMENT.md](AGREEMENT.md), root AGENTS.md.
Same real states and analytic kernel as the completed medium pilot. Only the
checkpoint changes scientifically. Pair/edge/node blocks are implementation
settings, preserving all neighborhoods and field couplings.

Protocol: `mace_polar_1l_analytic_multipole_vacuum_r2scan3c_pilot_v1`.
Campaign: `workspaces/mace_large_20260916/pilot_v2`.
Manifest SHA256: `fab9f9bcb96cbd06e1f268d595e18391c3c914ff3452ae52a582ed1ea0324fe6`.
Checkpoint SHA256: `9f65f8dc6ddaff1d631e299cb531376a7da5e68d1bef04f34a2d5073d5ef114b`.
Completed: job 1200525, all 12 calls, one A5000, 16 CPUs, 64,474 MiB requested
host RAM. Rotation/charge pass; partition fails. [REPORT.md](REPORT.md).
No executor remains active and no recomputation is needed to inspect results.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
MACE_PY="$PWD/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python"
LARGE_ROOT="$PWD/workspaces/mace_large_20260916/pilot_v2"
LARGE_MANIFEST="$LARGE_ROOT/manifest.json"
LARGE_RUNNER="$LARGE_ROOT/implementation/mace_hybrid.py"
"$MACE_PY" "$LARGE_RUNNER" dry-run --manifest "$LARGE_MANIFEST"
"$MACE_PY" "$LARGE_ROOT/implementation/mace_analytic_pilot.py" core-gate --manifest "$LARGE_MANIFEST"
"$MACE_PY" "$LARGE_RUNNER" collect --manifest "$LARGE_MANIFEST"
```

The existing executor preserves attempts and safely reuses verified successful
results. Eight core calls must pass charge and rotation checks before full calls.
The 2 kcal/mol partition diagnostic is reported separately. No compatible aquo
reference or calibration is available; missing scores remain null.

Preparation adds `--model-variant large --target-software PATH` to the existing
analytic-pilot preparer. The source collection is the completed analytic medium
`workspaces/mace_analytic_20260916/pilot_v1/collection_job_1200470.json`; target
software manifest is `workspaces/mace_large_20260916/software_v1/software_manifest.json`.
The isolated medium installation is reused without upgrades. Large's own
independent kernel checks are pinned by `kernel_test_receipt_v1.json`.

The legacy schema field `finite_reference` points to *analytic medium* in this
large manifest. Collected data explicitly label the comparison
`analytic_large_minus_analytic_medium`; no finite-displacement effect is mixed
into the model-size comparison. The medium and large results retain distinct
protocol/cache identities.

Terminal products are `collection_job_1200525.json`, workspace `REPORT.md`,
`job_1200525_accounting.json`, and per-task forces/densities/receipts. Preparation
v1 failed admission on a metadata-key mismatch before inference and is preserved;
v2 fixes only that check. Software records also preserve the initial DNS failure,
failed pre-load inspection and successful per-command resolution/download.
No shared resolver or installed dependency was changed.

Saved-output audit (no model inference), using a fresh output directory:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$MACE_PY" scripts/mace_output_audit.py \
  --medium-collection workspaces/mace_analytic_20260916/pilot_v1/collection_job_1200470.json \
  --large-collection "$LARGE_ROOT/collection_job_1200525.json" \
  --plan diagnostics/mace_large_20260916/OUTPUT_AUDIT_PLAN.md \
  --output workspaces/mace_large_20260916/output_audit_reproduction_v1
```

Existing `workspaces/mace_large_20260916/output_audit_v1/` contains its completed
result, mapped all-atom table and direct contrast gradient arrays. The command
refuses an existing output directory. It does not compute a hybrid gradient or
mechanical correction. Regression tests requiring local real artifacts skip
explicitly when those artifacts are absent.
