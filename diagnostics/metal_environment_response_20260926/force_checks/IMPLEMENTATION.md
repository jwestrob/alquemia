# Native force-check harness delivery — 27 September 2026

The finite 20-endpoint harness is implemented and passes six real-artifact tests.
No new molecular endpoints ran during implementation. Native derivative
qualification remains **pending execution**; the existing center results are
reused, not recomputed or presented as successful finite-difference checks.

Authoritative new-run agreement: `../FORCE_CHECK_PLAN.md`; original tolerances
come from `../PLAN.md`. This directory's PLAN.md records the mapped construction.
Root owns submission and the prepared production-of-research workspace. Nothing
changes the PQQ workflow, method defaults, old energies or source preparation.

## Operations

From repository root:

```bash
python scripts/metal_environment_force_checks.py prepare \
  --source-manifest workspaces/metal_environment_response_20260926/reference_scout_v1/manifest.json \
  --source-collection workspaces/metal_environment_response_20260926/reference_scout_v1/FINAL_COLLECTION_1219207.json \
  --agreement diagnostics/metal_environment_response_20260926/FORCE_CHECK_PLAN.md \
  --output workspaces/metal_environment_response_20260926/force_checks_v1 \
  --workers 20 --mpi-ranks 17
python workspaces/metal_environment_response_20260926/force_checks_v1/implementation/metal_environment_force_checks.py dry-run \
  --manifest workspaces/metal_environment_response_20260926/force_checks_v1/manifest.json
```

`execute` and `collect` accept the same manifest; `collect --output PATH` writes
an immutable result. Run execution only in the root-owned Slurm allocation,
with the appropriate PQQ dependency and resource layout. The example prepare
operation is implemented, not a claim that this named workspace already exists.
Do not rerun preparation over an existing directory.

The implementation snapshots native runner and renderer sources, including the
memory wrapper and unchanged base renderer. It retains exact template identity
and full source, state, method, ORCA, tolerance and implementation cache keys.
All source pin checks and generated coordinates are revalidated before execute.
Completed A centers are admitted only after existing runner receipt validation
against their original manifest, then actual output/gradient parsing.

## Geometry and interpretation

Thr159 hydrogen arc derivatives use unit physical tangents. Glu177 CA changes
its mapped cap through the exact source Jacobian and constant archived rounding
offset. Its removed MM charge is explicitly absent in the actual mapping, so
no MM atom moves in that mode. No cap is interpreted as an independent physical
motion. Both sides of the metal pair receive identical construction rules.

Collectors retain signed analytic and finite-difference derivatives at both
steps, step convergence, and repeat/rigid energy and separate QM/MM gradient
residuals. Any missing or failed cell blocks selected-direction qualification.
Even complete passing results would qualify only these derivatives of the
embedded electronic component, not missing hybrid terms or protein mechanics.

## Executed implementation checks

```bash
python -m unittest discover -s tests -p test_metal_environment_force_checks.py -v
```

Six passed, zero skips, 9.36 seconds for this run. Checks use real pinned maps,
archived native outputs/gradients/receipts, and a temporary complete manifest.
Tests cover exact repeat bytes, sole moved atoms, hydroxyl bond preservation,
analytic cap Jacobian finite differences, rigid distance invariance, corrupted
real mappings/collection status, and incomplete manifest rejection. Test
preparation writes temporary source-derived inputs only, then removes its own
temporary directory. No simulated successful molecular results or labels exist.
