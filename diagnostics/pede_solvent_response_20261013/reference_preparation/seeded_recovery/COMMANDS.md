# Conditional same-state Ca vacuum recovery — prepared only

**Closed without activation:** primary3372 subsequently converged and passed
collection. `READY.json` pins its actual successful receipt and collection.
The seeded contingency remains unrun. Do not launch its parser tests or molecular
calculation from the historical commands below; no recovery is now needed.

Scientific scope is root's `diagnostics/pede_solvent_response_20261013/SEEDED_RECOVERY_PLAN.md`.
At preparation, primary3372 was running; this preparation did not authorize a concurrent
retry or change its native iteration settings. The earlier description of125
as a total-iteration bound was incorrect; see the final electronic-limits note.
No molecule, parser test or Slurm job
was executed by the preparation agent.

Frozen workspace:
`workspaces/pede_solvent_response_20261013/reference_sample1_Ca_seeded_v1/`.

| Item | SHA-256 |
|---|---|
| manifest.json | `e4c700939bd4b7dae94de745fd5a76b33598dc6e3433fd20066f8b332c3d1848` |
| implementation/SEEDED.py | `b71ef49a43451d3acae2ac221703d3e0d8c68f116544aab4713c69b5ec10d26d` |
| source and staged initial.gbw | `1d718711e1aebe6720d8b4e5bbf0fab90cba9a930cc2aeb096f843a1c71421a3` |
| unchanged target XYZ | `2ac9319acde749d296132cccf71f45f87ede975fa3f4658dcdec2d68f509c07a` |

The one endpoint uses native ORCA6.1.1 r2SCAN-3c, DefGrid3, TightSCF, charge−3,
singlet838 electrons and the exact192-atom vacuum target. The14 caps and zero
explicit waters are retained. The target template differs only by an explicit
`%scf` block selecting `Guess MORead`, `MOInp "initial.gbw"` and `MaxIter125`.
`NoAutostart` remains present. The source GBW is copied and pinned; its successful
water energy is not reused as a vacuum energy. Original primary files, seven
completed endpoints and calibration are untouched.

The existing generic native runner and allocation-aware memory renderer are
copied unchanged. The current primary preflight is not relaxed or called on an
incompatible input. A narrow new preflight instead verifies the original source
and target manifests, source acceptance/receipt, atom order/coordinate bytes,
state, method/basis identity, source GBW completion timestamp and exact seed copy.
It reuses the original output collector plus the already tested TRAH parser fix.

## Small saved-artifact checks: one CPU

Root may supply healthy-node/output options. No memory/time override is needed.

```bash
SEED_WORK=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pede_solvent_response_20261013/reference_sample1_Ca_seeded_v1
SEED_SHA=e4c700939bd4b7dae94de745fd5a76b33598dc6e3433fd20066f8b332c3d1848
sbatch "$SEED_WORK/implementation/CHECK.sbatch" test \
  "$SEED_WORK/manifest.json" "$SEED_SHA" "$SEED_WORK/check_v1"
```

`test` includes preflight and six parser checks: a real executed MORead fixture,
rejection of an actual atomic guess, rejection of input-echo-only proof, a wrong
read filename, a missing basis-match message and model-potential fallback.
Negative fixtures exist only in memory. They establish parser behavior, not
Ca-vacuum convergence. Expected actual receipts are `check_v1/PREFLIGHT.json`,
`TESTS.json` and `EXECUTION.json`. Nothing currently asserts they passed.

For preflight alone, replace `test` with `dry-run` and use a fresh output folder.
For collection of an existing recovery, replace it with `collect` and another
fresh folder. All three actions require actual Slurm and make no molecular calls.

## Conditional native execution: twelve MPI ranks

After root's terminal-primary review and a passed test receipt, hash the actual
`check_v1/EXECUTION.json` and provide its digest in place of `ACTUAL_TEST_SHA`:

```bash
sbatch "$SEED_WORK/implementation/EXECUTE.sbatch" \
  "$SEED_WORK/manifest.json" "$SEED_SHA" "$SEED_WORK/execute_v1" \
  "$SEED_WORK/check_v1/EXECUTION.json" ACTUAL_TEST_SHA
```

The adapter currently requires a real **failed terminal native primary receipt**.
It refuses an absent receipt or a converged primary and acquires both workspace
locks, preventing a concurrent retry. A separately reviewed stop without such a
receipt is not silently treated as failure; that would need explicit preserved
termination evidence and a separately reviewed admission change. No automatic
retry, iteration extension or alternate spin/charge is implemented.

Execution retains actual allocation-derived MaxCore and native MPI binding.
Collection requires runtime `INITIAL GUESS: MOREAD`, the read-file line naming
`initial.gbw`, the matching-basis message, and absence of model-potential
fallback, all outside echoed input. It also checks printed vacuum/no-CPCM,
charge/spin, MaxIter125 and all existing native basis/electron/convergence/
population/wavefunction checks. It records any printed Aufbau reassignment;
same formal state does not mean fixed orbital occupations or ground-state proof.

Fresh outputs are `execute_v1/PREFLIGHT.json`, `EXECUTION_PREFLIGHT.json`,
`COLLECTION.json`, `EXECUTION.json`, plus native files under
`tasks/sample1_Ca_vacuum_seeded/`. The collection denominator is one conditional
recovery, explicitly separate from the original eight-cell primary. Root owns
submission, completion/health wakes and any later comparison.
