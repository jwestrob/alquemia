# Full100 triple transfer operations

All finite jobs, including the diagnostic and CPU preflight recoveries, are
complete and recorded under `run_v2` and its two `searches/shard_*` directories.
Do not submit duplicates. Final status is in CURRENT.md/REPORT.md and actual
receipts. These commands only inspect or
collect existing artifacts unless explicitly labeled as historical submissions.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
TRIPLE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$TRIPLE_PY" -m unittest discover -s tests -p test_union_triple_transfer.py -v
"$TRIPLE_PY" scripts/union_triple_transfer.py validate_origins \
 --manifest workspaces/union_triple_transfer_20260923/run_v2/origins/solvent/manifest.json
```

The fixed source selection, model/reference and all125 pair mappings are in
`run_v2/AUDIT.json`; finite110 search/cross/440 scalar requests are in
`run_v2/MOLECULAR_SCOPE.json`. Actual new origin forces are in
`run_v2/origins/COLLECTION.json`. All accepted pool reuses preserve their original
manifest/output/receipt pins. Preparation never creates substitute energies.

After both terminal pool collections exist, this computes the complete125-pair/
100-triple comparison without new molecular calls. Use a new output filename for
a replay; existing results are immutable.

```bash
"$TRIPLE_PY" scripts/union_triple_transfer_compare.py \
 --audit workspaces/union_triple_transfer_20260923/run_v2/AUDIT.json \
 --collections \
 workspaces/union_triple_transfer_20260923/run_v2/searches/shard_0/pool/collection_final.json \
 workspaces/union_triple_transfer_20260923/run_v2/searches/shard_1/pool/collection_final.json \
 --output workspaces/union_triple_transfer_20260923/COMPARISON_replay_v1.json
```

Origin/search/scalar wrappers, actual submission commands and their source hashes
are pinned in `SUBMISSION_origin_{gpu,cpu}.json` and
`searches/shard_*/SUBMISSION_{gpu,cpu}.json`. The immutable execution adapter is
`run_v2/execution_v1/union_triple_transfer_run.py` for searches; validator-only
host-copy recovery used `execution_v2/union_triple_transfer_run.py`. Each invocation supplies the
stage's implementation directory as PYTHONPATH. No input path requires source
editing. GPU search completion always leads to collection, including explicit
failed proposals; dependent CPU stages use afterany rather than hiding failures.

This is a frozen reference-source transfer experiment. It does not accept
arbitrary unknown proteins, invent canonical membership or change production.
