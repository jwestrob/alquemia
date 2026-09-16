# Matched archived baseline contrasts

**2026-09-16 — archive-only extraction and requested algebra; no new runs.**
The new global partition manifest pins the four old CPCM+MBIS outputs through
`baseline_output` and `source_pilot`. Their XYZ hashes match the new global
inputs exactly. Execution receipts, charge/multiplicity, native r2SCAN-3c,
CPCM(Water), DefGrid3, MBIS, ORCA 6.1.1 and native ECP/D4/gCP checks pass.

## Raw contrasts

`R = E_Ca − E_La`. Subtract Hartree endpoint energies first, then convert once
with the recorded factor **627.509474 kcal mol⁻¹ Hartree⁻¹**.

| Core | E_La, Hartree | E_Ca, Hartree | R, Hartree | R, kcal/mol |
|---|---:|---:|---:|---:|
| qm33 | −1755.178308184436 | −2401.246012601455 | −646.067704417019 | −405413.605367111069 |
| qm36 | −1983.765564490981 | −2629.829756462711 | −646.064191971730 | −405411.401274415315 |

**Partition difference, R(qm36) − R(qm33): +0.003512445289 Hartree =
+2.204092695754 kcal/mol.** Printed energy tokens and exact decimal arithmetic
are retained in the JSON; displayed precision is computational provenance,
not a claim of physical accuracy.

## Protocol and interpretation

- Endpoint protocol: `native_r2scan3c_cpcm_mbis_endpoints_v1`.
- qm33 source core: `pqq_vertical_swap_r2scan3c_native_cpcm_typed31_v2`.
- qm36 source core: `pqq_vertical_swap_r2scan3c_native_cpcm_typed31_qm36_v3`.

These are **matched baseline-derived pilot endpoints on the boundary-test
cores**, not the published canonical PQQ fixed-core calibration or holdout
values. No old APBS transfer term, aquo offset, calibrated zero or PQQ decision
band is included. The large common raw contrast has no classification meaning
by itself; the partition difference is the relevant comparison here. It does
not change the predeclared acceptance tolerance for the new global model.

Source pins, the four execution receipts and verification flags are in
[MATCHED_BASELINE_COMPARISON.json](MATCHED_BASELINE_COMPARISON.json), SHA256
`08088606734587ed2423b9e6ddaf2aad242b5818bd995f15c3bf77ccb3109d0b`.
Authority: `workspaces/global_electrostatic_20260916/partition_tasks_v1/manifest.json`
and its pinned `workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json`.
All four old outputs belong to completed job **1198934**.
