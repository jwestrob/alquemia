# Analytic MACE pilot operations

Standing pilot authorization is recorded in [AGREEMENT.md](AGREEMENT.md) and
root AGENTS.md. This is an opt-in research method; the baseline remains default.

Protocol: `mace_polar_1m_analytic_multipole_vacuum_r2scan3c_pilot_v1`.
Campaign: `workspaces/mace_analytic_20260916/pilot_v1`.
Manifest SHA256: `6b40455a9ca0cdf3249091cc12c722ff845005fdc23e72c827ce38841a494b72`.
Software: existing isolated `workspaces/mace_hybrid_20260916/software_v1`.
Exact launch command and resources are in the campaign's `launch.json`.
Final state and scientific findings belong in [REPORT.md](REPORT.md).

## Audit, gate, collect and report

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
MACE_PY="$PWD/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python"
ANALYTIC_ROOT="$PWD/workspaces/mace_analytic_20260916/pilot_v1"
ANALYTIC_MANIFEST="$ANALYTIC_ROOT/manifest.json"
ANALYTIC_RUNNER="$ANALYTIC_ROOT/implementation/mace_hybrid.py"
"$MACE_PY" "$ANALYTIC_RUNNER" dry-run --manifest "$ANALYTIC_MANIFEST"
"$MACE_PY" "$ANALYTIC_ROOT/implementation/mace_analytic_pilot.py" core-gate --manifest "$ANALYTIC_MANIFEST"
"$MACE_PY" "$ANALYTIC_RUNNER" collect --manifest "$ANALYTIC_MANIFEST"
```

The batch exit writes `collection_job_JOBID.json`. To export a new human report,
use `mace_hybrid.py report --collection PATH --output NEW_PATH`; writers refuse
to overwrite outputs. The completed campaign's full report is `REPORT.md` in
its workspace. No compatible aquo reference or calibrated class exists.

## Reproducible preparation and recovery

`scripts/mace_analytic_pilot.py prepare` takes explicit `--source-collection`,
`--rotation-manifest`, `--kernel-tests`, `--agreement` and `--output`. It copies
code into a new immutable campaign and preserves original input bytes. The
current inputs are the completed finite-displacement collection in
`workspaces/mace_hybrid_20260916/blocked_v6/collection_job_1200381.json`,
the original rotation manifest in
`workspaces/mace_rotation_20260916/audit_v1/original/manifest.json`, and
`workspaces/mace_analytic_20260916/kernel_test_receipt_v1.json`.

The kernel-test receipt pins the actual kernel, tests, log, checkpoint and
archived real densities. A stale or failed kernel receipt blocks admission.
Full-protein execution additionally requires all four primary/rotated core pairs
to pass the frozen numerical gates. The partition residual is reported
separately and is not concealed by the numerical gate.

The existing `run_pilot.sbatch` and MACE executor run this manifest. One A5000
share uses 16 CPUs, 64,474 MiB requested host RAM, `--partition=gpu`,
`--nodelist=node-128-512g-8gpu-1`, and
`--export=ALL,MACE_MIN_MEMORY_MIB=64474`. The cluster seven-day QOS limit remains;
there is no project time/CPU stopping budget. No active executor should be
started twice. Only interrupted/failed tasks are retried; successful receipts
are verified and reused. Read the terminal collection before launching recovery.
