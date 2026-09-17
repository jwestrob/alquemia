# Explicit field + MACE short component: operations

Run from the repository root. Inputs, source states, software, checkpoint,
parameters and eight tasks are pinned in the manifest; no production change.
Preparation directory creation is exclusive. Choose a new output directory for
an intentional preparation replay. Never overwrite scientific outputs.

```bash
ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FIELD_WORK=$ALCH/workspaces/mace_omol_20260917
FIELD_MANIFEST=$FIELD_WORK/explicit_field_short_v1/manifest.json
cd "$ALCH"

# Audit / validate exactly the pinned preparation; no scientific computation.
"$FIELD_PY" "$FIELD_WORK/explicit_field_short_v1/implementation/mace_explicit_field_short.py" dry-run --manifest "$FIELD_MANIFEST"

# Preparation replay: use a new, nonexistent output path.
"$FIELD_PY" scripts/mace_explicit_field_short.py prepare \
  --gb "$FIELD_WORK/full_boundary_GB_report_v1/result.json" \
  --charges "$FIELD_WORK/normalized_charge_v1/manifest.json" \
  --reuse "$FIELD_WORK/field_short_reuse_audit_v1.json" \
  --short-reference "$ALCH/workspaces/mace_short_engine_20260916/pilot_v1/collection_job_1200736.json" \
  --agreement "$ALCH/diagnostics/mace_omol_20260917/EXPLICIT_FIELD_SHORT_PLAN.md" \
  --output "$FIELD_WORK/explicit_field_short_replay_v1"

# Already submitted as 1200980 and 1200981. Do not duplicate live jobs.
# Accepted results are reused on restart; failed utility attempts require
# explicit --retry-failed inside an existing allocation.
sbatch --output="$FIELD_WORK/explicit_field_short_v1/potential_%j.log" \
  diagnostics/mace_omol_20260917/run_explicit_potential.sbatch "$FIELD_MANIFEST"
sbatch --partition=gpu --nodelist=node-128-512g-8gpu-1 --gres=gpu:1 \
  --cpus-per-task=16 --mem=64474M --export=ALL,MACE_MIN_MEMORY_MIB=64474 \
  --job-name=alquemia-field-short \
  --output="$FIELD_WORK/explicit_field_short_v1/short_%j.log" \
  diagnostics/mace_hybrid_20260916/run_pilot.sbatch \
  "$ALCH/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" \
  "$FIELD_MANIFEST" native

# Collect, compare and report from real accepted receipts. With no --output,
# this is read-only. An incomplete input never becomes a zero correction.
"$FIELD_PY" "$FIELD_WORK/explicit_field_short_v1/implementation/mace_explicit_field_short.py" collect --manifest "$FIELD_MANIFEST"
```

No new DFT, charge fitting, whole-protein inference, solvent solve or training.
The eight short readouts include their own analytic gradients; these are not
combined hybrid forces. Exact density potential uses immutable copied native
wavefunctions. Point-charge coupling is a representation diagnostic only.
See the plan for energy accounting, frozen gates and evidence limitations.
