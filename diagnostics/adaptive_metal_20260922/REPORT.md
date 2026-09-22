# Joint metal response is justified for a bounded test; preparation is ready

The completed angular searches leave appreciable native-MACE forces on the
frozen metal in all eight endpoints. A separate seven-coordinate proposal
adapter is prepared on the same four consumed contexts. **No optimization,
MACE, solvent or DFT calculation was launched.** No prediction changed.

| Source | Ca metal gradient norm | La metal gradient norm |
|---|---:|---:|
| 1H4I | 21.461652 | 40.269650 |
| 4MAE | 27.054198 | 38.961880 |
| PLM8344 sample0 | 16.500426 | 31.451218 |
| PLM07ab sample0 | 13.332442 | 24.154770 |

Units are kcal/mol/Å; each value is evaluated at that endpoint's actual final
angular candidate. These are **vacuum-MACE gradients**, not a matched-geometry
Ca−La difference or solvent-composite gradient. Seven selected angular spaces
have residual loads below 0.0007 kcal/mol/Å; PLM8344 La retains its final physical
boundary and 2.459946 selected angular load. The residual metal forces motivate
the additional coordinate block, but do not establish a useful correction.

## Prepared implementation

`scripts/adaptive_metal_proposals.py` reuses the existing Kinematics, native warm
MACE worker, source/cap/overlap checks and constrained SLSQP settings. The three
metal coordinates precede the exact same original-q0-selected four angles.
Both metals have identical selected IDs. The physical displacement constraint
applies to every source heavy atom and constrains the metal **norm** to 0.8 Å;
three component bounds alone would permit an excessive diagonal displacement.

Objective derivatives carry explicit per-coordinate units: eV/Å for the metal,
eV/rad for the angles. Force-to-gradient sign and cap chain rules remain explicit.
There is no all-radian field for the mixed vector. Each endpoint minimizes its
own native energy from its original q0; no favorable starting point, relative-score
objective or new chemical state is introduced. Water/PQQ/scaffold coordinates are
fixed. Failed searches retain unavailable candidates and null solvent/score fields.

The older broad 0.20 Å/0.20 rad metal-plus-all-donor experiment already worsened
native DFT separation. This is a controlled addition to the newer force-selected
four-angle subspace with its existing 0.8 Å final domain. It may also fail. A
shared candidate-pool comparison is still required to establish any utility.

## Real checks and current numerical limitation

Seven tests pass in **16.875 s, zero skips**: all eight real maps/source pairs,
common selectors, explicit mixed units, the seven-column constraint Jacobian,
archived Cartesian-force projection with correct sign/units/caps, the physical
metal sphere, fixed nonselected atoms and unavailable pre-execution outputs.
Geometry-only differences of a frozen archived force covector test the mapping;
they are not new finite-difference MACE energies. No scientific backend was mocked.

The current SLSQP interface can abort on an oversized trial that triggers the
unchanged severe-overlap guard. The actual MMOL1770 continuation failure is
diagnosed in [SEARCH_FAILURE.md](SEARCH_FAILURE.md). **This behavior is preserved
in the prepared joint version**, pending root review of a separately declared
numerical change. No guard was loosened and no energy was assigned to a clash.

## Artifacts

- Manifest: `workspaces/adaptive_metal_20260922/prepared_v1/manifest.json`
- SHA-256: `3276a4927093d08e6d780d9b951389d01d5aadf7b93c05261a710c7a275e7324`
- Exact consumed force receipts/vectors: `prepared_v1/RESIDUALS.json`
- Complete unrun eight-endpoint denominator: `prepared_v1/collection_unrun_v1.json`
- Frozen plan: [PLAN.md](PLAN.md); commands: [COMMANDS.md](COMMANDS.md)
- Test log: [TESTS_v1.txt](TESTS_v1.txt)

New scientific allocation cost: **zero**. Root reviews before any calls. Baseline,
production/default, original angular pilot, continuation and pool engine unchanged.
