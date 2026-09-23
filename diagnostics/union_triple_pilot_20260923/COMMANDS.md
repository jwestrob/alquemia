# Approved threefold-preparation pilot operations

All paths explicit, original chemistry/defaults preserved. Live finite jobs are
recorded in CURRENT.md; do not duplicate them. CPU automatic collection follows
terminal GPU state, including retained failures.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
TRIPLE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$TRIPLE_PY" -m unittest discover -s tests -p test_union_triple_pilot.py -v
"$TRIPLE_PY" workspaces/union_triple_pilot_20260923/proposals_v1/implementation/union_triple_adaptive.py validate \
 --manifest workspaces/union_triple_pilot_20260923/proposals_v1/manifest.json
```

Actual origin collection is `origins_v1/COLLECTION.json`. The proposed and then
actual shared pool is `pool_v1/manifest.json`; terminal collection is
`pool_v1/collection_final.json`. Final29-source rows and the distinct canonical25
reference are in `REFERENCE_v1.json`. These are under
`workspaces/union_triple_pilot_20260923/`.

Read-only scientific collection replay after the jobs finish, writing a new file:

```bash
"$TRIPLE_PY" workspaces/union_triple_pilot_20260923/pool_v1/implementation/union_triple_adaptive.py collect_pool \
 --manifest workspaces/union_triple_pilot_20260923/pool_v1/manifest.json \
 --output workspaces/union_triple_pilot_20260923/pool_v1/collection_replay.json
"$TRIPLE_PY" workspaces/union_triple_pilot_20260923/pool_v1/implementation/union_triple_adaptive.py reference \
 --collection workspaces/union_triple_pilot_20260923/pool_v1/collection_replay.json \
 --output workspaces/union_triple_pilot_20260923/REFERENCE_replay.json
```

The precise preparation/submission commands and all implementation pins are in
`origins_v1/SUBMISSION_{gpu,cpu}.json` and `SUBMISSION_adaptive_{gpu,cpu}.json`.
Those commands are execution receipts, not instructions to rerun completed work.
No source editing is needed for explicit request/preparation paths. The finite
adapter intentionally rejects a different source count/optimizer profile. Full
primary transfer is not included in this execution interface.
