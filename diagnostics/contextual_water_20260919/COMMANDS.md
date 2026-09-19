# Explicit contextual-water operations

Run from repository root. Outputs are write-once; the commands below use fresh
named replay destinations. Existing `replay_v3` and `score_replay_v1` are complete.

```bash
CW_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$CW_PY" scripts/contextual_water_prepare.py dry-run \
  --manifest workspaces/contextual_water_20260919/replay_v3/manifest.json
"$CW_PY" scripts/contextual_water_score.py dry-run \
  --manifest workspaces/contextual_water_20260919/score_replay_v1/manifest.json
"$CW_PY" -m unittest discover -s tests -p test_contextual_water.py -v
```

Create a request for explicitly supplied prepared cases (this concrete example
uses the two real alpha structures). Omit `--context-group` for independent cases.
The input interface itself needs no case-name changes in source code.

```bash
"$CW_PY" scripts/contextual_water_prepare.py request \
  --preparation workspaces/mace_global_benchmark_20260916/prepared_v1/ALPHA_1F6S/preparation.json \
  --preparation workspaces/mace_global_benchmark_20260916/prepared_v1/ALPHA_6IP9/preparation.json \
  --context-group ALPHA_1F6S=alpha_same_indexed_construct \
  --context-group ALPHA_6IP9=alpha_same_indexed_construct \
  --policy-source diagnostics/hydration_network_20260918/CONFIG.json \
  --proposal-source-manifest workspaces/hydration_network_20260918/mace_optimization_v1/manifest.json \
  --scorer-source-manifest workspaces/hydration_scanner_20260919/mace_v2/manifest.json \
  --agreement diagnostics/contextual_water_20260919/AGREEMENT.md \
  --output workspaces/contextual_water_20260919/explicit_request_v1.json
"$CW_PY" scripts/contextual_water_prepare.py prepare \
  --request workspaces/contextual_water_20260919/explicit_request_v1.json \
  --output workspaces/contextual_water_20260919/explicit_prepare_v1
"$CW_PY" scripts/contextual_water_prepare.py dry-run \
  --manifest workspaces/contextual_water_20260919/explicit_prepare_v1/manifest.json
```

The simple request above intentionally has no cached-result pins: it renders four
pending proposal tasks. To reuse the demonstrated actual computations, use the
recorded `diagnostics/contextual_water_20260919/REPLAY_REQUEST.json` as the request
instead. Reuse pins are checked against coordinates, state, method and receipts.

Execution must occur inside a suitable existing GPU allocation. The same
execute/collect interface handles partial jobs and skips accepted completed tasks:

```bash
"$CW_PY" scripts/contextual_water_prepare.py execute \
  --manifest workspaces/contextual_water_20260919/explicit_prepare_v1/manifest.json
"$CW_PY" scripts/contextual_water_prepare.py collect \
  --manifest workspaces/contextual_water_20260919/explicit_prepare_v1/manifest.json \
  --output workspaces/contextual_water_20260919/explicit_prepare_v1/collection.json
"$CW_PY" scripts/contextual_water_score.py prepare \
  --collection workspaces/contextual_water_20260919/explicit_prepare_v1/collection.json \
  --output workspaces/contextual_water_20260919/explicit_score_v1
"$CW_PY" scripts/contextual_water_score.py dry-run \
  --manifest workspaces/contextual_water_20260919/explicit_score_v1/manifest.json
"$CW_PY" scripts/contextual_water_score.py execute \
  --manifest workspaces/contextual_water_20260919/explicit_score_v1/manifest.json
"$CW_PY" scripts/contextual_water_score.py collect \
  --manifest workspaces/contextual_water_20260919/explicit_score_v1/manifest.json \
  --output workspaces/contextual_water_20260919/explicit_score_v1/collection.json
"$CW_PY" scripts/contextual_water_score.py report \
  --collection workspaces/contextual_water_20260919/explicit_score_v1/collection.json \
  --output workspaces/contextual_water_20260919/explicit_score_v1/report
```

Do not rerun the already-completed alpha demonstration merely to exercise this
example. This task executed only cached replay and zero-work execute operations.
The CPU interpreter prepares/collects source chemistry; the runner uses the pinned
MACE interpreter for scientific workers. No DFT calls are required by this component.
