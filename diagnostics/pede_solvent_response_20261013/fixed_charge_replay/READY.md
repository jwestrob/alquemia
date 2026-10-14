# Native fixed-charge ALPB replay — ready, unexecuted

Frozen manifest:
`workspaces/pede_solvent_response_20261013/fixed_charge_replay_v1/manifest.json`

SHA256 `2469b38c0d227df8806aa48763608c50672eb1660a14532679c65185c214b5ca`.

No compiler, molecular backend, replay or scientific test ran on login.
`PREPARE.py` only joined existing metadata and pinned files. Root submits and
arms the actual completion wake. The previous seeded DFT fallback stays closed.

## One command

```bash
REPLAY_WORK=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pede_solvent_response_20261013/fixed_charge_replay_v1
REPLAY_SHA=2469b38c0d227df8806aa48763608c50672eb1660a14532679c65185c214b5ca
sbatch --nodelist=node-64-768g-5 "$REPLAY_WORK/implementation/RUN.sbatch" \
  build-test-execute "$REPLAY_WORK/manifest.json" "$REPLAY_SHA" \
  "$REPLAY_WORK/execution_v1"
```

Standard shared allocation,12 task CPUs, no memory/time override. Compilation
uses at most12 workers. Each of the eight component calls uses one numerical
thread. No SCC, DFT, gradients, optimization or environment installation occurs.
`execution_v1/FINAL_COLLECTION.json` is the terminal watcher result, including
failure receipts. Existing output directories are refused.

## Native implementation and precision

The unchanged xTB6.7.1 factory loads its internal GFN2-water parameters. The
wrapper calls `newBornModel`, `TBorn%update`, and `TBorn%getEnergyParts` directly.
It requires the actual internal parameter label. Named settings reproduce the
archived P16/ALPB, GBOBC, gsolv,298.15K,230-grid, zero-salt setup; radii, solvent
dielectric and empirical coefficients come from the native factory. Electronic
smearing300K is already represented in the saved SCC endpoints and is not a
new solvent thermodynamic temperature.

Sources are immutable archives of xTB6.7.1 commit26b28010e805f7d1aeeef39813feb473e69cc4be,
mctc-lib0.4.2 commit51139755818af265803e600059cf3c44fe14995e, and
test-drive0.4.0 commit9c3401e30dbd2da1add77aaa252a4c6928fe39a1. The latter is a
build dependency; its tests are not run. The installed native package also uses
mctc-lib0.4.2. Optional tblite/CPCM-X backends are disabled. Actual compiler,
OpenBLAS, build commands, resulting binary and loaded libraries are retained.

Actual native charge artifacts store eight decimals. No binary restart or
higher precision is invented. `qsh` is a zero-length argument because the exact
pinned component routine does not reference shell charges. With rounding radius
`d=5e-9e`, the frozen bounds are:

- Gelec: `d*sum(abs(B_e*q)) + d²/2*sum(abs(B_e))`, with native electrostatic
  `B_e=bornMat−diag(2*hbw)`.
- Ghb: `2*d*sum(abs(hbw*q)) + d²*sum(abs(hbw))`.
- SASA/shift: no charge-rounding contribution.

Each receives the prior1e-9Eh native component-closure allowance. These matrix
operations bound uncertainty only; all reported energies come from native
`getEnergyParts`. A tolerance exceeding the existing GFN backend parity1e-6Eh
is explicitly uninformative and fails admission. No thresholds depend on replay
outcomes. Geometry, atom order, charge file, accepted result and source state
are cross-checked in Slurm.

## Finite gate and accounting

Four water-charge evaluations must reproduce every archived Gelec/Ghb/Gsasa/
Gshift component within its predeclared bound. Only then run the four matching
vacuum-charge evaluations. Same-case cavity and native matrix hashes must be
identical across charge distributions. The exact archived intrinsic remainder
`I` is retained: `U=G(qvac)+[G(qwater)−G(qvac)+I]` within the diagonal rounding
budget. The collection includes raw endpoint values and fixed double contrasts
`sample1_Ca−sample1_La−sample3_Ca+sample3_La` for both terms and their sum,
with one archived Hartree-to-kcal conversion. This is a response diagnostic,
not a new solvent correction or classifier.

If root deliberately separates stages, the same wrapper accepts `build`, then
`test --build-receipt BUILD.json`, then
`execute --build-receipt BUILD.json --qualification QUALIFICATION.json` with
distinct fresh output folders. Reuse the completed diagonal receipts; do not
repeat the combined command after success. No automatic retry is implemented.
