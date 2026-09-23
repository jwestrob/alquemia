# Consistent-context research operations

This branch is opt-in. The production scanner and all historical references stay
unchanged. Existing completed calculations need no resubmission. All paths below
are explicit; each command writes new output and refuses to overwrite a record.

```bash
UNION_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
UNION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Read the frozen reference and completed full transfer

```bash
"$UNION_PY" - "$UNION_ROOT" <<'PY'
import json, pathlib, sys
root = pathlib.Path(sys.argv[1])/'workspaces/consistent_context_20260922'
reference = json.loads((root/'calibration28_v1/REFERENCE_v2.json').read_text())
print('Frozen bands:', reference['bands'], 'gap:', reference['gap_model_kcal_mol'])
result = json.loads((root/'primary225_v1/comparison_v1.json').read_text())
print(json.dumps(result['counts']['single_sources']['all'], indent=2))
PY
```

## Reproduce geometry preparation only

The actual original source manifest and fixed membership selection are reused;
this operation performs no molecular energy calculation.

```bash
"$UNION_PY" "$UNION_ROOT/scripts/consistent_context.py" prepare \
  --preparation "$UNION_ROOT/workspaces/accommodation_goal_20260920/folds_v1/preparation_reconciled_v1.json" \
  --agreement "$UNION_ROOT/diagnostics/consistent_context_20260922/PLAN.md" \
  --output "$UNION_ROOT/workspaces/consistent_context_20260922/preparation_replay_v1" \
  --workers 4
```

## Validate the finite full-transfer manifest

```bash
"$UNION_PY" "$UNION_ROOT/workspaces/consistent_context_20260922/primary225_v1/implementation/consistent_context.py" validate-solvent \
  --manifest "$UNION_ROOT/workspaces/consistent_context_20260922/primary225_v1/solvent/manifest.json"
```

The prepared `READY.json` records fresh and reused task counts. The per-phase
submission JSON files record exact commands and file hashes. In the full stage,
`submission_mace_layout_v2.json` and `submission_collect_layout_v2.json` are the
executed replacements of cancelled, never-started requests; the scientific
manifests are unchanged. `run_mace.sbatch` calls the
existing warm-MACE executor; `run_solvent.sbatch` calls the existing ORCA runner.
The dependent collector performs no molecular calls and preserves failed cells.
All jobs are now complete. See REPORT.md, RESULT_v1.json and COST_v1.json.

## Re-collect and compare saved results without scientific execution

Run after the declared molecular jobs finish. Automatic collection produces the
original `collection_final_v1.json` and `comparison_v1.json`; the commands below
write a distinct replay. They preserve the frozen canonical-only reference.

```bash
"$UNION_PY" "$UNION_ROOT/scripts/consistent_context_compare.py" collect \
  --stage "$UNION_ROOT/workspaces/consistent_context_20260922/primary225_v1" \
  --output "$UNION_ROOT/workspaces/consistent_context_20260922/primary225_v1/collection_replay_v1.json"

"$UNION_PY" "$UNION_ROOT/scripts/consistent_context_compare.py" transfer-compare \
  --scope primary225 \
  --collection "$UNION_ROOT/workspaces/consistent_context_20260922/primary225_v1/collection_replay_v1.json" \
  --reference-path "$UNION_ROOT/workspaces/consistent_context_20260922/calibration28_v1/REFERENCE_v2.json" \
  --baseline "$UNION_ROOT/workspaces/nikasha_recovery_20260922/proposal_comparison.json" \
  --sources "$UNION_ROOT/diagnostics/accommodation_controls_20260920/PQQ_ALL250_SOURCES.json" \
  --output "$UNION_ROOT/workspaces/consistent_context_20260922/primary225_v1/comparison_replay_v1.json"
```

## Actual-fixture checks

```bash
cd "$UNION_ROOT"
"$UNION_PY" -m unittest discover -s tests -p 'test_consistent_context*.py' -v
```

No test substitutes fabricated energies for unavailable scientific executables.
Malformed-input checks explicitly corrupt copies of actual archived artifacts.
