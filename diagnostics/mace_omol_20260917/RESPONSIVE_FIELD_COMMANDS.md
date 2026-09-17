# Responsive field candidate: explicit operations

All commands run from the repository root. Scientific outputs are immutable;
use a new output directory for a replay. Do not duplicate live allocations.
The physical plan is RESPONSIVE_FIELD_PLAN.md; native accounting is
RESPONSIVE_FIELD_ACCOUNTING.md. Production baseline/default remains unchanged.

```bash
ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FIELD_WORK=$ALCH/workspaces/mace_omol_20260917
cd "$ALCH"

# Eight native quantum endpoints: completed job1200983. No rerun is needed.
# Audit / dry-run of their immutable scientific inputs:
"$FIELD_PY" scripts/mace_responsive_field.py dry-run \
  --manifest "$FIELD_WORK/responsive_quantum_v2/manifest.json"

# Recollect from saved outputs; the preserved collection parser recognizes the
# exact conditional ORCA warning while still rejecting numerical gradients.
"$FIELD_PY" "$FIELD_WORK/responsive_quantum_collection_v2/implementation/mace_responsive_field.py" collect \
  --manifest "$FIELD_WORK/responsive_quantum_v2/manifest.json"

# Full replay preparation, if independently needed, without execution:
"$FIELD_PY" scripts/mace_responsive_field.py prepare \
  --source "$FIELD_WORK/explicit_field_short_result_v1.json" \
  --agreement "$ALCH/diagnostics/mace_omol_20260917/RESPONSIVE_FIELD_PLAN.md" \
  --accounting "$ALCH/diagnostics/mace_omol_20260917/RESPONSIVE_FIELD_ACCOUNTING.md" \
  --output "$FIELD_WORK/responsive_quantum_replay_v1"
# Execution uses the existing64CPU /4x16rank wrapper:
# sbatch diagnostics/mace_omol_20260917/run_responsive_quantum.sbatch \
#   "$FIELD_WORK/responsive_quantum_replay_v1/manifest.json"
```

Quantum result `responsive_quantum_result_v2.json` is accepted; V1 retained the
initial parser rejection. V1preparation and V2execution have identical physical
inputs; the first preparation was superseded to preserve an old parser message.
The collector repair required no new SCF. Point-charge inputs, source mappings,
charge/multiplicity and method remain fixed.

Chargepreflight1200984ran no fitting/potential calls. Worker diagnostic1200985
proved only1.11e-16cap-weight roundoff; V3 preserves recorded weights and allows
1e-12arithmetic replay. Same scientific tolerances and inputs. Utilities1200986
completed8fits+8potentials; allfit/projection/coupling gates pass. Solver1200989
contains44newcalls and4matched environment-only reuses; inspect status first.

```bash
# Actual charge preparation audit / restart (only if no live utility job).
"$FIELD_PY" "$FIELD_WORK/responsive_charges_v3/implementation/mace_responsive_charges.py" dry-run \
  --manifest "$FIELD_WORK/responsive_charges_v3/manifest.json"
# sbatch diagnostics/mace_omol_20260917/run_normalized_charges.sbatch \
#   "$FIELD_WORK/responsive_charges_v3/manifest.json"

# Native-fit replay from accepted embedded quantum outputs, no execution.
"$FIELD_PY" scripts/mace_responsive_charges.py prepare \
  --quantum "$FIELD_WORK/responsive_quantum_result_v2.json" \
  --output "$FIELD_WORK/responsive_charges_replay_v1"

# A report writes exclusively to a new directory. Existing accepted report is
# responsive_charge_report_v2/result.json; V1 is the preserved incomplete stage.
"$FIELD_PY" "$FIELD_WORK/responsive_charges_v3/implementation/mace_responsive_charges.py" report \
  --manifest "$FIELD_WORK/responsive_charges_v3/manifest.json" \
  --output "$FIELD_WORK/responsive_charge_report_replay_v1"

# Prepare a new solver replay from the actual accepted new density.
"$FIELD_PY" scripts/mace_responsive_solvent.py prepare \
  --charges "$FIELD_WORK/responsive_charge_report_v2/result.json" \
  --software "$ALCH/workspaces/mace_gb_20260916/pilot_v2/software.json" \
  --output "$FIELD_WORK/responsive_GB_replay_v1"

# Audit original prepared solver inputs.
FIELD_GB=$FIELD_WORK/responsive_GB_v1/manifest.json
"$ALCH/workspaces/mace_gb_20260916/software_v2/venv/bin/python" \
  "$FIELD_WORK/responsive_GB_v1/implementation/mace_hybrid.py" dry-run --manifest "$FIELD_GB"

# Execution already submitted as1200989. Never duplicate a live allocation.
# sbatch --partition=gpu --nodelist=node-128-512g-8gpu-1 --gres=gpu:1 \
#   --cpus-per-task=16 --mem=64474M --export=ALL,MACE_MIN_MEMORY_MIB=64474 \
#   diagnostics/mace_hybrid_20260916/run_pilot.sbatch \
#   "$ALCH/workspaces/mace_gb_20260916/software_v2/venv/bin/python" "$FIELD_GB" native

# Read-only collect + comparison. Missing/failed terms stay unavailable.
"$FIELD_PY" "$FIELD_WORK/responsive_GB_v1/implementation/mace_responsive_solvent.py" collect --manifest "$FIELD_GB"
```
