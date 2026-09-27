# Finite classical component checks — executable, not yet evaluated

The harness implements the frozen `../../COMPONENT_CHECK_PLAN.md` on the exact
real ledger. **No molecular energy/force calls ran during implementation.**
Four real-artifact construction/mapping/accounting tests pass in 18.998 s.
OpenMM Context creation is patched to fail in those tests, preventing an
accidental molecular evaluation. One initial construction test caught an
optional `kind` field on protein particles; fixed before any production snapshot.

## Executable accounting

`scripts/metal_environment_component_checks.py` provides prepare, dry-run,
execute and collect. Preparation creates two immutable OpenMM XML Systems,
16 coordinate arrays shared by the two metals, and a 32-task manifest. It
snapshots its own implementation and the common artifact helper.

The four force groups are:

0. Native ff19SB retained bonds, angles, torsions and CMAP, reconstructed from
   pinned parent force objects and selected original term indices. All original
   CMAP tables are copied. No constrained dynamics or integrator steps occur.
1. MM–MM LJ through native NonbondedForce; real QM epsilon is zero here.
   Every original protein exception replaces the generic pair; QM or cross
   exception contributions are explicitly zero in this group.
2. MM–MM Coulomb through native NonbondedForce using qtilde. Real QM charges
   are zero only in this *MM-only* force object; this is not a model of QM
   charges. Native exception pair products use the ledger's corrected values.
3. Real QM–MM LJ through OpenMM CustomNonbondedForce with one interaction
   group, Lorentz–Berthelot mixing, and all protein exceptions excluded.
   A CustomBondForce supplies the exact 67 cross exception replacements.
   All PQQ/metal noncovalent pairs use the declared parameter family.

Both native and custom nonbonded objects use NoCutoff. No dispersion-tail or
reaction-field term, C4 induction, cap FF particle, classical QM–QM interaction,
or classical QM–MM Coulomb is present. Custom forces use OpenMM's derivative
engine, not a new numerical engine. Physical positions convert Å→nm once.
OpenMM quantities convert energy to kcal/mol and forces to kcal/mol/Å once.

## Manifest and qualification

Per metal: A, B, exact repeat A, joint rigid A, then three physical directions
at ±0.001/±0.0005 Å. Directions are the exact-source Thr159 hydroxyl arc,
Glu177 CA displacement perpendicular to CB→CA, and metal translation along
metal→Glu177 OE1 fixed at A. Metal and boundary motions move real atoms only;
no cap is an independent classical coordinate. Full 9,141-atom forces are
retained for every component/configuration in compressed NumPy artifacts.

Collection reports 200 checks: 32 repeat/rigid component residuals, 72 FD/step
convergence residuals, and 96 paired-metal MM-component equality checks across
all 16 geometries. Paired equality uses the declared repeat-energy and rigid-force
limits. Every FD has its analytical projection and raw central-difference value.
No failed cell is silently replaced or omitted; complete qualification requires
all 32 cells and all checks. The force is negated for the energy derivative.

Archive successful receipts unchanged. Their manifest hash binds state/metal,
geometry, ledger, model, implementation, numerical limits and software version.
Required pins are checked before execution; each receipt records its own task,
model, Slurm job, wall time, component timings, task CPU time, and process peak
RSS. Peak RSS is explicitly a process-lifetime maximum, not per-task memory.
A previously successful scout is reused. Failed receipts are preserved and need
an explicitly new attempt directory; there is no automatic retry policy.

Only Slurm execution creates Contexts. Reference platform is explicit;
ProcessPoolExecutor bounds independent single-thread work to the requested
worker count and actual allocation. Contexts remain warm within each worker
when it handles repeated cells of a metal. Root first measures Ca_A on one
shared CPU, then decides whether to run the remaining 31; the script does not
submit or automatically expand the scout.

This qualifies classical component derivatives only. The complete dry-system
hybrid is not qualified, earlier electronic rotational failures remain, and
solvent is absent. No optional quantum-gradient combination or optimization is
implemented. Solvent-consistent whole-protein relaxation remains unavailable.

## Commands

From the repository root:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/metal_environment_component_checks.py prepare \
  --ledger workspaces/metal_environment_response_20260926/hybrid_preparation_v1/ledger_v1/LEDGER.json \
  --agreement diagnostics/metal_environment_response_20260926/COMPONENT_CHECK_PLAN.md \
  --output workspaces/metal_environment_response_20260926/component_checks_v1
```

After preparation, use the frozen script:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  workspaces/metal_environment_response_20260926/component_checks_v1/implementation/metal_environment_component_checks.py \
  dry-run --manifest workspaces/metal_environment_response_20260926/component_checks_v1/manifest.json
```

**Inside root's one-CPU Slurm allocation only**, replace `dry-run` with:

```text
execute --manifest workspaces/metal_environment_response_20260926/component_checks_v1/manifest.json --workers 1 --task-id Ca_A
```

The separately scheduled remaining matrix uses `execute ... --workers 32`
and reuses completed Ca_A. Collect with the same frozen script:

```text
collect --manifest workspaces/metal_environment_response_20260926/component_checks_v1/manifest.json --output workspaces/metal_environment_response_20260926/component_checks_v1/FINAL_COLLECTION.json
```

Implementation tests (no energy calls):

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  -m unittest discover -s tests -p test_metal_environment_component_checks.py -v
```

Use a fresh output filename for each partial/final collection. The harness never
rewrites a scientific collection or an existing preparation directory.
