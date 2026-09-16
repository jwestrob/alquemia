# Approved global electrostatic experiment: operations

These commands refer to the frozen approved artifacts. They do not add cases,
change parameters or authorize Stage 2 from partial results. Current quantum
and ESP work is complete. Surface jobs 1199956 and 1199959 are running; do not
submit duplicate jobs while either corresponding slice is active.

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
squeue -j 1199956,1199959
sacct -j 1199949,1199952,1199956,1199959 --parsable2 \
  --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,TotalCPU,MaxRSS,NodeList
```

The initial job owns only `1h4i_qm33_La/primary` and
`1h4i_qm33_Ca/primary`. The second slice owns the other 23 tasks. Every solver
uses one CPU; independent tasks fill the requested shared CPU allocation.

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

sbatch --ntasks=23 --mem=48G \
  --output="$FIELD_WORK/surface_recovery_%j.out" \
  --error="$FIELD_WORK/surface_recovery_%j.err" \
  diagnostics/global_electrostatic_20260916/run_surfaces.sbatch \
  "$FIELD_WORK/surfaces_remaining_execution_v1/campaign_manifest.json"
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
