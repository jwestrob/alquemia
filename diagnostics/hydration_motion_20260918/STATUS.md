# Completed — 2026-09-18

All four states pass the local radial response criteria. DFT1201953 completed8
endpoints; MACE1201958 completed16calls. No outstanding motion job or next-stage
submission. Read REPORT.md, RESULT_v2.json, RECEIPTS.json and COMMANDS.md.
Baseline/PQQ unchanged;56distinct real-fixture tests pass. Entropy/occupancy
remain unavailable pending coupled physical-coordinate and basin validation.
Total69,008allocated core-s/137GPU-s includes startup failure. Vault updated.

Native orientation comparator1201824 has completed all four1F6S optimizations.
They converge normally but fail the frozen-coordinate tolerance (drift2.26e-5 to
7.46e-4Å); one also exceeds the rigid-water tolerance. Actual final energies and
gradient traces are retained in orientation_1f6s_v1/collected_opt_v2.json, without
promoting them to qualified exact-geometry minima. No rerun or tolerance change.
Comparator1201825 for6IP9 remains live; preserve it. These are separate older
runs and do not affect the successful exact-coordinate motion checks.

## Retained launch checkpoint

# Water-motion checkpoint — 2026-09-18

Approved by Jacob: “proceed.” Scope is AGREEMENT.md. New native DFT job1201953
ran eight radial-displacement EnGrad endpoints (four16-rank workers,64CPU).
MACE job1201958 ran16 identical/finer displaced energy-force calls on oneGPU,
16CPU,64474MiB. Four actual center energies/gradients are reused. No production
or PQQ changes and no occupancy/entropy correction is enabled by this test.

Workspaces: DFT/preparation in hydration_motion_20260918/pilot_v2;
MACE retry in hydration_motion_20260918/mace_retry_v1, referencing the exact
same preparation and XYZ files. Initial pilot_v1 failed a path-containment dry
run before submission; pilot_v2 copies identical DFT XYZ bytes under its manifest.
MACE1201954 failed before inference because exact equality of a recomputed unit
vector differed by2.22e-16 across NumPy versions. Retry checks the physical vector
at1e-14 roundoff tolerance and validates displaced coordinates exactly against
the pinned original vector. Scientific acceptance tolerances are unchanged.
No calculated result was discarded or recomputed. Preserve failed records.

Native orientation comparators1201824/1201825 are separate prior work; preserve
and collect when complete. Do not resubmit them or old completed occupancy jobs.
