# Exact four-cell continuation operations

Job1211126 executes the fixed two-pass diagnostic. Do not resubmit completed
cells. The submission, wrapper, finite stage1 manifest, exact seeds and agreement
are pinned in `SUBMISSION.json`. Stage2 is prepared from actual confirmed
stage1 seed pairs; failures remain in the four-cell denominator.

Read-only validation of the executed snapshot:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/precision_geometry_continuation_20260923/run_v2/stage1/implementation/precision_geometry_continuation.py \
validate --manifest workspaces/precision_geometry_continuation_20260923/run_v2/stage1/manifest.json
```

Once both collections exist, optional comparison replay to a new output:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/precision_geometry_continuation_20260923/run_v2/stage2/implementation/precision_geometry_continuation.py \
compare --stage1 workspaces/precision_geometry_continuation_20260923/run_v2/stage1/collection.json \
--stage2 workspaces/precision_geometry_continuation_20260923/run_v2/stage2/collection.json \
--output workspaces/precision_geometry_continuation_20260923/run_v2/COMPARISON_replay_v1.json
```

Preparation v1 was never submitted. Its only subsequent code correction lets
an entirely failed first stage produce a stage2 manifest with zero executable
tasks and four missing rows. The executed v2 snapshot preserves that failure
accounting; recipe, seeds, source coordinates and scientific scope are unchanged.
