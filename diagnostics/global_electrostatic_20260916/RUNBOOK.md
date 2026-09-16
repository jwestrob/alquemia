# Approved global electrostatic experiment: operations

These commands refer to the frozen approved artifacts. They do not add cases,
change parameters or authorize Stage 2 from partial results. Current quantum
and ESP work and initial surface pair1199956 are complete. Recovery array1199964
indices0–8 and group1199974 are active. Original batch1199959 was stopped after a documented persistent
low-clock observation; its completed isolated controls and partial attempts
remain. Do not submit duplicate work while corresponding tasks are active.

## Paths and executables

Run this setup in the shell used for the commands below:

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
FIELD_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
FIELD_WORK="$FIELD_ROOT/workspaces/global_electrostatic_20260916"
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

- ORCA: `/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca`.
- ESP utility: the same directory's `orca_vpot`.
- Pinned solver copies: `workspaces/global_electrostatic_20260916/software/tabi_native_v1/`;
  `software_manifest.json` records the exact TABI/NanoShaper executable paths
  and hashes. The adapter reads this manifest; no PATH-selected substitute.
- Quantum inputs: `partition_tasks_v1/manifest.json`.
- Authoritative solver schedule: `surfaces_v1/campaign_manifest.json`.
- Remaining-task execution slice: `surfaces_remaining_execution_v1/campaign_manifest.json`.
- Current recovery: array indices0–8 in `surface_recovery_array_v1/array_manifest.json`,
  plus `surface_recovery_group_v1/campaign_manifest.json` for the other12 tasks.

Output names below are explicit and new. Writers refuse an existing output;
retain earlier versions and use a new recorded version for another collection.

## Audit and dry-run

Read [AGREEMENT.md](AGREEMENT.md), [NUMERICS.md](NUMERICS.md),
[BOUNDARY_AUDIT.md](BOUNDARY_AUDIT.md), and [ACCURACY_INPUTS.md](ACCURACY_INPUTS.md).
These checks launch no QM or solvent calculation:

```bash
"$FIELD_PY" scripts/global_electrostatic.py dry-run \
  --manifest "$FIELD_WORK/partition_tasks_v1/manifest.json" \
  --output "$FIELD_WORK/partition_tasks_v1/dry_run_operator_v1.json"

"$FIELD_PY" workspaces/global_electrostatic_20260916/accuracy_input_audit_v1/audit_inputs.py \
  --output "$FIELD_WORK/accuracy_input_audit_v1/operator_recheck_v1.json"

"$FIELD_PY" -m unittest discover -s tests -p 'test_global_electrostatic_accounting.py' -v
```

## Collect completed quantum endpoints

The verified collection already exists at
`partition_tasks_v1/collection_verified_v1.json`. Recollection reads existing
outputs and execution receipts; it does not rerun ORCA:

```bash
"$FIELD_PY" scripts/global_electrostatic.py collect \
  --manifest "$FIELD_WORK/partition_tasks_v1/manifest.json" \
  --output "$FIELD_WORK/partition_tasks_v1/collection_operator_v1.json"
```

Keep new MBIS charges paired with their own gas-phase ESP receipts in `esp_v1`.
The old CPCM charge/ESP records cannot satisfy these tasks.

## Monitor and recover the existing surface work

```bash
squeue -j 1199964,1199974
sacct -j 1199949,1199952,1199956,1199959,1199964,1199974 --parsable2 \
  --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,TotalCPU,MaxRSS,NodeList
```

The initial job owns only `1h4i_qm33_La/primary` and
`1h4i_qm33_Ca/primary`. Recovery array1199964 indices0–8 own nine full-protein
tasks; its explicit slices are in `surface_recovery_array_v1/array_manifest.json`.
The other twelve pending array elements were cancelled without execution and
replaced by group1199974. Its wrapper verifies that the allocation spans at
least twelve distinct physical cores before launching the twelve workers.
The two isolated controls from1199959 are reused. Every solver uses one thread;
actual affinity/topology are retained for the recovery. No frequency, affinity,
queue-priority or node-administration setting was mutated.

**Recovery only, after the corresponding job has ended:** the same exact
commands below can resume those approved slices after a technical interruption.
The executor validates pins, takes campaign/task locks, reuses verified success
and preserves attempts. Do not remove locks or edit cached artifacts. Inspect
failed receipts before retrying; a physical/numerical failure does not authorize
a parameter change or repeated scientific rescue.

```bash
sbatch --ntasks=2 --mem=16G \
  --output="$FIELD_WORK/surface_recovery_%j.out" \
  --error="$FIELD_WORK/surface_recovery_%j.err" \
  diagnostics/global_electrostatic_20260916/run_surfaces.sbatch \
  "$FIELD_WORK/surfaces_v1/campaign_manifest.json" \
  1h4i_qm33_La/primary 1h4i_qm33_Ca/primary

sbatch --array=0-8 --exclude=node-48-256g-13 \
  --output="$FIELD_WORK/surface_array_recovery_%A_%a.out" \
  --error="$FIELD_WORK/surface_array_recovery_%A_%a.err" \
  diagnostics/global_electrostatic_20260916/run_surface_array.sbatch \
  "$FIELD_WORK/surface_recovery_array_v1/array_manifest.json"

sbatch --exclude=node-48-256g-13 \
  --output="$FIELD_WORK/surface_group_recovery_%j.out" \
  --error="$FIELD_WORK/surface_group_recovery_%j.err" \
  diagnostics/global_electrostatic_20260916/run_surface_group.sbatch \
  "$FIELD_WORK/surface_recovery_group_v1/campaign_manifest.json"
```

No CPU-time or project wall-time stopping budget is imposed. The native GMRES
iteration/convergence rule is part of the frozen solver. Slurm partition limits
and job receipts remain recorded; do not interpret a pending estimate as progress.

## Collect, assess and report the physical gate

Always collect from the **master 25-task manifest**, not the 23-task execution
slice. Incomplete collection remains explicitly incomplete. These commands
read existing artifacts, perform the approved accounting/checks, and submit
nothing:

```bash
"$FIELD_PY" scripts/global_electrostatic_assess.py collect \
  --campaign "$FIELD_WORK/surfaces_v1/campaign_manifest.json" \
  --output "$FIELD_WORK/surfaces_v1/collection_final_v1.json"

"$FIELD_PY" scripts/global_electrostatic_assess.py assess \
  --collection "$FIELD_WORK/surfaces_v1/collection_final_v1.json" \
  --output "$FIELD_WORK/surfaces_v1/assessment_final_v1.json"

"$FIELD_PY" scripts/global_electrostatic_assess.py report \
  --assessment "$FIELD_WORK/surfaces_v1/assessment_final_v1.json" \
  --output "$FIELD_WORK/surfaces_v1/physical_report_final_v1.md"
```

Retain the numerical report beside receipts. Update the concise tracked
[REPORT.md](REPORT.md) with the completed gate and measured costs, including
failed attempts and preparation. No automatic Stage 2 submission follows any
command above. Root must review the complete physical gate and the exact five
environment preparations; unsupported cases remain explicit.
