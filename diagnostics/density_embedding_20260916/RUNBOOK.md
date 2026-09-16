# Density and interaction diagnostics: operator commands

Run from the Alquemia repository root. These are opt-in research operations;
the production scorer, released PQQ bands and old experiments remain unchanged.
The [density report](REPORT.md), [charge-fit result](CHARGE_FIT_RESULT.md) and
[interaction scope](INTERACTION_PLAN.md) describe distinct experiments.

```bash
DIAGNOSTIC_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
DIAGNOSTIC_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
DIAGNOSTIC_WORK="$DIAGNOSTIC_ROOT/workspaces/density_embedding_20260916"
```

## Collect existing calculations without rerunning them

Each output must be new; writers refuse overwrites. Failed or incomplete native
calculations remain unavailable. The collectors check pinned inputs, execution
receipts, electronic states, coordinate mappings and energy accounting.

```bash
"$DIAGNOSTIC_PY" scripts/density_embedding.py collect-potentials \
  --manifest "$DIAGNOSTIC_WORK/potential_retry_v1/potential_manifest.json" \
  --output "$DIAGNOSTIC_WORK/potentials_operator_v1.json"
"$DIAGNOSTIC_PY" scripts/density_embedding.py collect-embedded \
  --manifest "$DIAGNOSTIC_WORK/prepared_v1/manifest.json" \
  --potentials "$DIAGNOSTIC_WORK/potentials_operator_v1.json" \
  --output "$DIAGNOSTIC_WORK/embedded_operator_v1.json"
"$DIAGNOSTIC_PY" scripts/density_embedding.py collect-chargefit \
  --manifest "$DIAGNOSTIC_WORK/chelpg_v1/manifest.json" \
  --output "$DIAGNOSTIC_WORK/chargefit_operator_v1.json"
"$DIAGNOSTIC_PY" scripts/interaction_decomposition.py collect \
  --manifest "$DIAGNOSTIC_WORK/eda_order_retry_v1/manifest.json" \
  --output "$DIAGNOSTIC_WORK/interaction_operator_v1.json"
"$DIAGNOSTIC_PY" scripts/interaction_decomposition.py report \
  --result "$DIAGNOSTIC_WORK/interaction_operator_v1.json" \
  --output "$DIAGNOSTIC_WORK/interaction_operator_v1.md"
```

The EDA collector distinguishes native convergence, component closure and
agreement with the original fragment references. Its native generated ghost
basis centers are not physical atoms. A converged EDA calculation does not
automatically pass the predeclared reference checks or diagnose the original
partition discrepancy. No affinity score or classification is produced.

## Reproduce preparation and dry-run

This creates the exact declared two-task Asp303 comparison, including source
atom mappings. It launches no calculation. Both fragments and adduct retain
native r2SCAN-3c/DefGrid3; the geometry declaration precedes `%Frag` as required
by ORCA. The first failed input-order attempt remains in `eda_v1/`.

```bash
"$DIAGNOSTIC_PY" scripts/interaction_decomposition.py prepare \
  --states-manifest "$DIAGNOSTIC_ROOT/workspaces/global_electrostatic_20260916/states_v1/states_manifest.json" \
  --endpoint-manifest "$DIAGNOSTIC_ROOT/workspaces/global_electrostatic_20260916/partition_tasks_v1/manifest.json" \
  --plan "$DIAGNOSTIC_ROOT/diagnostics/density_embedding_20260916/INTERACTION_PLAN.md" \
  --agreement "$DIAGNOSTIC_ROOT/diagnostics/density_embedding_20260916/AGREEMENT.md" \
  --output "$DIAGNOSTIC_WORK/eda_operator_prepared_v1"
"$DIAGNOSTIC_PY" scripts/affordable_workflow.py dry-run \
  --manifest "$DIAGNOSTIC_WORK/eda_operator_prepared_v1/manifest.json" \
  --output "$DIAGNOSTIC_WORK/eda_operator_prepared_v1/dry_run.json"
```

The existing `run_interaction.sbatch` accepts one absolute manifest path and
uses two 16-rank workers on 32 allocated CPUs, with one thread per rank. Original
submission commands are recorded in each attempt's `submission.json`. Native
EDA performs additional fragment SCFs internally; count those, not just the two
top-level tasks. It has no project CPU/time stopping budget. Completed outputs
are reused only with matching provenance; partial attempts require an explicit
new retry directory.

## Software verification

```bash
"$DIAGNOSTIC_PY" -m unittest discover -s tests -p test_density_embedding.py -v
"$DIAGNOSTIC_PY" -m unittest discover -s tests -p test_interaction_decomposition.py -v
```

These tests parse actual saved scientific outputs or prepare real frozen
structures. Deliberately corrupted copies test rejection. Missing real fixtures
are explicit skips; no synthetic successful scientific output substitutes for
an unavailable executable or calculation.
