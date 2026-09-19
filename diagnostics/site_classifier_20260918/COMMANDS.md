# Archived-feature classifier commands

No command here launches ORCA, MACE inference, a GPU job or model fine-tuning.
The scientific settings are fixed in AGREEMENT.md; changing them requires a
new experiment version. Every output must be a new path.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
CLASSIFIER_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
CLASSIFIER_WORK="$PWD/workspaces/site_classifier_20260918"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

# Reproduce features from pinned, already computed sources.
"$CLASSIFIER_PY" scripts/site_classifier.py inventory \
  --root "$PWD" \
  --agreement "$PWD/diagnostics/site_classifier_20260918/AGREEMENT.md" \
  --output "$CLASSIFIER_WORK/inventory_replay_v1"
"$CLASSIFIER_PY" scripts/site_classifier.py features \
  --inventory "$CLASSIFIER_WORK/inventory_replay_v1/inventory.json" \
  --output "$CLASSIFIER_WORK/features_replay_v1"

# Exact same declared fits, using the completed original feature table.
"$CLASSIFIER_PY" scripts/site_classifier.py evaluate \
  --features "$CLASSIFIER_WORK/features_v2/features.json" \
  --output "$CLASSIFIER_WORK/evaluation_replay_v1"

# Export existing fitted coefficients with held-out evaluation status; no refit.
"$CLASSIFIER_PY" scripts/site_classifier.py export \
  --result "$CLASSIFIER_WORK/evaluation_v1/result.json" \
  --output "$CLASSIFIER_WORK/models_replay_v1"

# Score a real available case with the exported research model.
"$CLASSIFIER_PY" scripts/site_classifier.py score \
  --model "$CLASSIFIER_WORK/models_v1/PQQ_functional_class__DFT_structure.json" \
  --features "$CLASSIFIER_WORK/features_v2/features.json" \
  --case 1H4I --output "$CLASSIFIER_WORK/score_1h4i_replay_v1.json"

"$CLASSIFIER_PY" tests/test_site_classifier.py
```

The three model names are `DFT`, `DFT_structure`, and `DFT_structure_MACE`.
Targets are `PQQ_functional_class` and `direct_site_affinity_direction`.
`score` rejects mixing targets or DFT/MACE protocols. Unlabelled rows can be
scored but cannot be counted as correct/incorrect. The direct-affinity exports
are explicitly unevaluable with the current held-out groups; the MACE-augmented
PQQ model is marked as a held-out regression. No artifact is the production default.

For another already computed case, supply its pinned preparation/endpoints in
the inventory schema and run `features`; no source-code path edit is needed.
Adding a new biological evaluation or energy calculation is separate from
replaying this completed experiment.
