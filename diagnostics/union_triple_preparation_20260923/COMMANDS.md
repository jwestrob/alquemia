# Three-fold membership preparation operations

These operations reconstruct source graphs and audit archived receipts only.
They launch no molecular calculations. Use new output paths for new records.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
TRIPLE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$TRIPLE_PY" -m unittest discover -s tests -p 'test_union_triple*.py' -v
```

## Existing completed artifacts

- Primary frozen selection and all100 triple rows:
  `workspaces/union_triple_preparation_20260923/run_v3/SELECTION.json`
- All125 deduplicated primary contexts plus the separate stress probe:
  `workspaces/union_triple_preparation_20260923/run_v3/PREPARATION.json`
- Exact native/solvent receipt mapping and explicit missing cells:
  `workspaces/union_triple_preparation_20260923/run_v3/REUSE.json`
- Canonical deterministic triple choices and skipped incomplete choices:
  `workspaces/union_triple_preparation_20260923/canonical28_v1/SELECTION.json`
- Canonical25 plus3 unchanged crystal preparations and reuse:
  `workspaces/union_triple_preparation_20260923/canonical28_v1/{PREPARATION,REUSE}.json`

## Explicit rerunnable preparation commands

```bash
"$TRIPLE_PY" scripts/union_triple_preparation.py select \
 --inventory workspaces/union_subset_inventory_20260923/inventory_v1.json \
 --output workspaces/union_triple_preparation_20260923/review_v1
"$TRIPLE_PY" scripts/union_triple_preparation.py prepare \
 --manifest workspaces/union_triple_preparation_20260923/review_v1/SELECTION.json \
 --workers 4
"$TRIPLE_PY" scripts/union_triple_preparation.py reuse \
 --preparation workspaces/union_triple_preparation_20260923/review_v1/PREPARATION.json \
 --tenfold-collection workspaces/consistent_context_20260922/primary225_v1/collection_final_v1.json \
 --local-inventory workspaces/accommodation_goal_20260920/folds_v1/INVENTORY.json \
 --local-comparison workspaces/accommodation_goal_20260920/folds_v1/comparison_v1.json \
 --output workspaces/union_triple_preparation_20260923/review_v1/REUSE.json
"$TRIPLE_PY" scripts/union_triple_reference.py \
 --primary workspaces/union_triple_preparation_20260923/review_v1/PREPARATION.json \
 --crystals workspaces/consistent_context_20260922/crystal_controls_v1/CRYSTALS.json \
 --output workspaces/union_triple_preparation_20260923/review_canonical_v1 \
 --workers 4
"$TRIPLE_PY" scripts/union_triple_preparation.py reuse \
 --preparation workspaces/union_triple_preparation_20260923/review_canonical_v1/PREPARATION.json \
 --tenfold-collection workspaces/consistent_context_20260922/calibration28_v1/collection_final_v2.json \
 --local-inventory workspaces/accommodation_goal_20260920/folds_v1/INVENTORY.json \
 --local-comparison workspaces/accommodation_goal_20260920/folds_v1/comparison_v1.json \
 --output workspaces/union_triple_preparation_20260923/review_canonical_v1/REUSE.json
```

`prepare` preserves all unsupported members and checks source atom identities,
bonds, H parents, caps, charge, cofactor/water inventory and paired coordinates.
The archived scalar/force mapping checks actual checkpoint and native recipes;
old eight-rank receipts keep their original numerical provenance. Candidate-pool
reuse is explicitly unaudited. This is not a scoring or submission CLI, and no
new classification or calibration is emitted.
