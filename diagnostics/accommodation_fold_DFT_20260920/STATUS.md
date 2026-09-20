# Fold challenge DFT baseline: running

2026-09-20. Frozen [plan](PLAN.md); no production/default changes.

**416 new native DFT endpoints submitted**, with25 verified canonical pairs
reused. All250 source geometries remain in the ledger:233 supported and17
unsupported; the primary pool is225 noncanonical sources,208 supported.
The existing failed MMOL hydrogen preparation remains unsupported here; its
separate numerical repair was not introduced into this experiment.

| Shard | New endpoints | Slurm job | Requested resources |
|---|---:|---|---|
| 0 | 104 | 1203976 | 64CPU,600GiB; four16-rank endpoints |
| 1 | 104 | 1203977 | 64CPU,600GiB; four16-rank endpoints |
| 2 | 104 | 1203978 | 64CPU,600GiB; four16-rank endpoints |
| 3 | 104 | 1203979 | 64CPU,600GiB; four16-rank endpoints |

All four started on separate64CPU nodes;16 actual endpoint outputs began
with ORCA6.1.1 and no stderr startup errors. The original canonical pairs
reproduce25/25 released bands.

First verified partial collection (`startup_collection.json`): **11/416
endpoints complete,405 pending, no recorded failures**. Five new pairs are
complete, all from A0A3F2YLY8: Ca-conditioned samples0/1/3 are inconclusive;
Ca4 and La3 are La-supported. This is an early within-protein robustness
observation, not a final panel comparison. The other220 primary-source scores
remain unavailable in that snapshot. Completed endpoint wall sum4,299.387923 s
at16 ranks (68,790.206768 endpoint-rank-seconds); partial-completion timing
does not establish final throughput.

## Checks actually run

- Native method/state matching for every supported source; exact copied
  inputs/XYZ and Ca/La coordinate equality; canonical coordinate agreement
  to the predeclared1e-12 Å reconciliation tolerance.
-50 real archived outputs, energies, execution receipts, runtime inputs and
  released bands verified for reuse; same executable hash verified locally.
-416-task dry-run: four nonoverlapping104-endpoint manifests.
- Prelaunch collection correctly reports416 pending endpoints,25 canonical
  results and225 unavailable primary scores; no zeros or fallback scores.
- Four focused real-fixture tests passed, including a deliberately corrupted
  real method input. Log [TESTS_FIXED.txt](TESTS_FIXED.txt). The initial test
  invocation had a test-harness attribute-name collision, repaired before
  submission; its traceback remains in [TESTS.txt](TESTS.txt).
- Join to the real MACE `native_comparison_v2.json` passes exact frozen
  source/preparation and label/group identity checks. No new model fitting.

## Artifacts and collection

Run: `workspaces/accommodation_fold_DFT_20260920/run_v1/`.
Master `manifest.json` SHA256:
`2d15d781feb0f32f5f0e3c2539cad05496a16f5fff79f1b21f4d4d06593245b5`.
Four shard manifests, submission receipts and a preserved five-file executor
implementation are under this run. Every scientific input is an unchanged
copy of the emitted fold core file.

Each wrapper collects even after a nonzero executor exit, preserving that
exit code. `collection_job_<id>.json/.md` reports all presently available
results and failures. Once all416 task receipts are terminal, the collector
also writes immutable `final_collection.json/.md`. Missing receipts remain
pending, including after a killed allocation; they are never interpreted as
computed failures or successful scientific endpoints. Any such interruption
needs explicit reconciliation from Slurm/runner logs before restart.

## Runnable operations (repository root)

```bash
# Verify the finite, already-submitted manifest; executes no quantum work.
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/accommodation_fold_DFT_20260920/run_v1/implementation/accommodation_fold_dft.py dry-run --manifest workspaces/accommodation_fold_DFT_20260920/run_v1/manifest.json

# Collect a new immutable snapshot after the jobs finish or for partial status.
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/accommodation_fold_DFT_20260920/run_v1/implementation/accommodation_fold_dft.py collect --manifest workspaces/accommodation_fold_DFT_20260920/run_v1/manifest.json --output workspaces/accommodation_fold_DFT_20260920/run_v1/collection_after_jobs.json --final-if-complete

# Match final DFT results to the frozen native-MACE comparison on identical sources.
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/accommodation_fold_DFT_20260920/run_v1/implementation/accommodation_fold_dft.py compare --dft workspaces/accommodation_fold_DFT_20260920/run_v1/final_collection.json --mace workspaces/accommodation_goal_20260920/folds_v1/native_comparison_v2.json --output workspaces/accommodation_fold_DFT_20260920/run_v1/final_native_MACE_comparison.json
```

The comparison accepts an explicit later MACE composite collection through
the same `--mace` option, provided its preparation/source records match.
Current MACE native/core, native/context and optional composite decisions
remain separate; no missing composite correction is replaced with zero.
Collection automatically writes the DFT Markdown count table. Counts cover
singlefolds, completeLa4/Ca5 medians and the equal mean of both arm medians.
Structural repeats are not independent biological observations, ensemble
descriptors are not thermal populations, and original single-source bands
are transferred without refitting.

## Cost evidence

All50 archived compatible endpoints:17,822.282511 s summed wall,
349.3619455 s median and267,334.237665 endpoint-rank-seconds, with15 ranks
on node-344-8t-1. New jobs use16 ranks on64CPU nodes; these historical costs
are not a matched-hardware throughput comparison. New endpoint receipts
retain actual wall/parallelism/allocation. No new GPU calls or custom budget.
