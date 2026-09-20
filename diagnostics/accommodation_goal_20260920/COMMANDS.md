# Fold robustness experiment: actual commands

Run from the repository root on biotite. Do not resubmit completed work. The actual
manifest, implementation and output hashes are saved below workspaces/accommodation_goal_20260920/folds_v1.
All paths are explicit; fresh commands require a fresh output directory.

```bash
export FOLD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

Actual preparation used prepare_folds.sbatch, job1203744:16one-thread processes,
250sources, unchanged seeded chemistry. Native execution used run_fold_mace.sbatch,
job1203769 with mace_v3/manifest.json: one H200,32CPUs,200000MiB. The warm calculator
retains analytic forces and resets its ASE cache for each new physical endpoint.
The first2canonical endpoints serve as the existing0.01kcal execution check;
both match exactly. Remaining98canonical endpoints are reused with explicit
state/model/source-coordinate checks. No GPU tensor batching approximation.

Actual solvent submissions and exact executable commands are in
[FOLD_SUBMISSIONS.json](FOLD_SUBMISSIONS.json). Four disjoint source-case manifests
retain the existing8MPI-rank/concurrent8-task runner per64CPU node. Only their
execution grouping changed; input bytes, scientific keys and1864endpoint scope
are preserved. Each job writes its collection_JOBID.json automatically.

Read the actual frozen source plan and source data with:

```bash
cat diagnostics/accommodation_goal_20260920/FOLD_ROBUSTNESS_PLAN.md
cat diagnostics/accommodation_goal_20260920/FOLD_SUBMISSIONS.json
```

Preparation/adapters and restartable native endpoint operation are exposed by:

```bash
"$FOLD_PY" scripts/accommodation_folds.py --help
"$FOLD_PY" scripts/accommodation_folds.py prepare --help
"$FOLD_PY" scripts/accommodation_folds.py prepare-mace --help
"$FOLD_PY" scripts/accommodation_folds.py run-mace --help
"$FOLD_PY" scripts/accommodation_folds.py split-solvent --help
"$FOLD_PY" -m unittest discover -s tests -p test_accommodation_folds.py -v
```

Comparison/reporting is in scripts/accommodation_folds_compare.py. Its actual
commands and final results will be recorded after the four solvent collections
exist. The original25replay cases stay separate from225otherstructural samples.
No threshold is refitted, missing scores remain missing, and no ensemble is
silently accepted after dropping one of its declared members.

The fast production source-backed interface is documented separately in
[released commands](../pqq_fast_release_20260920/COMMANDS.md).
