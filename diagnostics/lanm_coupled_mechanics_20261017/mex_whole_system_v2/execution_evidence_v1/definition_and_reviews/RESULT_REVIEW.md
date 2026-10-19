# Actual preflight 14789 review

**PASS for the declared single before-La, multiplicity-10 H200 scout.** This is
an input and saved-solvent admission result; no whole-system MACE energy or force
has been evaluated.

The completed result reports `preflight_passed` on `node-64-768g-1`, with zero
model loads, native forwards, molecular calls and optimizer steps. All eight
energy/force result slots remain explicitly unavailable/not attempted.

## Inspected evidence

- Result: `workspaces/lanm_coupled_mechanics_20261017/mex_whole_preflight_v2/execution/RESULT.json`,
  SHA256 `7a33f499c6eb5fe696c65432902324e0b0886f64398c45be0d315dc508f0fdaa`.
- Manifest: same directory, `MANIFEST.json`,
  SHA256 `acdaa4623a85a6e36fa7edc23d43a94fd94bc51ff9465aea69834a41f1b8d8d0`.
- Executed source: `dd611ca03982e4f2548a88dabf9477e33aba7592fc2b41d832e60d13edd2add7`;
  executed config: `ee9f5173063751a62f7c66d962549b1d2568c611f8053aff54d9280eb43da598`.

The manifest preserves the exact job12970 start/terminal and job12981 paired
cells: two geometries × La/Dy × two aggregate-spin hypotheses. Every task has
2,049 atoms, charge +2 and target index2045. The target and all three Nd
spectators retain their source IDs; all167 waters are retained, with no caps or
MM/GFN interaction terms added. Reported electron counts are7805 for La and7814
for Dy. The declared multiplicities are La10/Dy15 and La4/Dy9; preflight's
unchanged parity and checkpoint-category checks passed.

All four saved ALPB admissions include2,049 physical gradient rows. Their
receipt, XYZ, charge, gradient, native execution, START and stdout pins were
inspected. Exact component-gradient replay passed during the real preflight.
The before-La/before-Dy/after-La/after-Dy receipt SHA prefixes remain
`7c812b2b`, `bcb67b34`, `49073ec0`, `3f5ce4a1`.

The actual three negative controls rejected a changed spectator atom, changed
whole charge and inherited local-core spin. Result and manifest scientific
bindings match.

## Limits and accounting

Reported application elapsed time:2.688187s; measured process CPU:2.346512s;
peak parent RSS:61,304KiB; allocation:one CPU. These are application records,
not independently inferred scheduler utilization.

The first native scout must still establish finite full forces, actual graph
state and H200 resource feasibility. Remaining seven calls require its original
600-second forward and90%-of-device-memory gates. The aggregate-spin categories
do not identify local coupling or establish physical state or model accuracy.

Original failed preflight14788 remains preserved. This review read saved
artifacts and performed pin/metadata checks only; no numerical rerun or
scientific import was performed.
