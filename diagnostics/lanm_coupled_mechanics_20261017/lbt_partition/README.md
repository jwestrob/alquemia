# Executable LBT local/MM preparation

Source-only delivery. Root owns submission. AST parsing and shell syntax checks
passed; no scientific imports, numerical preparation or energy calls were run
on the login node.

From the repository root, the prepared wrapper accepts three absolute paths:

```bash
sbatch diagnostics/lanm_coupled_mechanics_20261017/lbt_partition/RUN.sbatch \
  "$PWD/diagnostics/lanm_coupled_mechanics_20261017/lbt_partition/PREPARE.py" \
  "$PWD/diagnostics/lanm_coupled_mechanics_20261017/lbt_partition/CONFIG.json" \
  "$PWD/workspaces/lanm_coupled_mechanics_20261017/lbt_partition_v1/execution"
```

The output directory must be fresh. The wrapper requests one CPU and no memory
or time override. It uses the pinned existing `lanm_qmmm` Python/OpenMM8.5.1.
No `Context`, MACE, xTB, ORCA, solvent call or optimizer is created.

## Outputs and scope

`RESULT.json` retains four source/region slots and eight declared metal states,
including any unavailable preparation. `REGION_SELECTION.json` records the
same source-union selection for both structures, before any new energy.

Each prepared source/region contains:

- `PARTITION.json`, `MAPPING.json`, `PHYSICAL.json`, `core_mapping.json` and
  `environment_atoms.json`: exact physical membership and current cap/field map.
- `La_Q.xyz`, `Dy_Q.xyz`, `environment.pc`: paired local geometry and shifted field.
- `LEDGER.json`, `particles.json`, `exceptions.json`, `bonded_terms.json`,
  `BONDED_TERM_DECISIONS.json`: the complete native interaction accounting.
- `La_MM.xml`, `Dy_MM.xml`: actual native systems containing only retained MM
  terms; these have not been energy/force evaluated.
- `GLOBAL_SOLVENT_CHARGES.json`: separate unshifted whole-physical charges.
- `PREPARATION.json`: state, provenance, geometry-map derivative checks and the
  distinction between preparation and unqualified molecular forces.

The region policy is in [REGION_POLICY.md](REGION_POLICY.md). Empty MM is exactly
zero, with no empty interaction group. A baseline that already covers the whole
peptide has no distinct expanded region; that slot stays unavailable.

## Reusable coordinate contract

Within a Slurm allocation, load this exact module and call `initialize(cfg)`
with its frozen configuration, then read an admitted `PARTITION.json`.

```python
mapping = update_geometry(partition, physical_positions_A)  # shape (256, 3)
physical_gradient = scatter(mapping, qm_gradient, pc_gradient)
```

Both input gradients must use the same units; output preserves those units and
the positive-energy-gradient convention. No force-sign or unit conversion is
implicit. Every cap updates from both true physical anchors and the field follows
its physical source IDs. All 256 physical coordinates remain represented.

Use the serialized MM systems with the current physical coordinates. The future
complete evaluator must qualify its own electronic/MM/solvent derivative and
bind actual electronic executables/parameters; these preparation task rows are
not permission to use defaults. This model replaces whole-peptide MACE with a
local electronic term plus the explicit environment ledger. Never sum both.
