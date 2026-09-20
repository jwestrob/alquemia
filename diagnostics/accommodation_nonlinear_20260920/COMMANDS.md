# Bounded nonlinear donor pilot

Job1204162 completed the immutable pilot_v2 manifest. Do not resubmit or edit it.
All commands below run from the repository root with the existing CPU environment.

Read the final report and compact actual result:

```bash
cat diagnostics/accommodation_nonlinear_20260920/REPORT.md
```

The full final collection is
`workspaces/accommodation_nonlinear_20260920/collection_final_v1.json`.
All eight candidates are retained; zero complete pairs pass minimum qualification.
The separate native-candidate execution gate remains closed.

Validate the pinned finite experiment (no molecular calls):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 scripts/accommodation_nonlinear.py validate \
 --manifest workspaces/accommodation_nonlinear_20260920/pilot_v2/manifest.json
```

Replay the completed endpoints into a new immutable collection, without launching
any calculation. This explicit output path is reserved for a future read-only replay:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 scripts/accommodation_nonlinear.py collect \
 --manifest workspaces/accommodation_nonlinear_20260920/pilot_v2/manifest.json \
 --output workspaces/accommodation_nonlinear_20260920/collection_replay_v1.json
```

The batch wrote final `result_1204162.json`, outer execution
receipt and warm-GPU summary under pilot_v2. Each endpoint retains its optimizer
result and all successful/failed coordinate evaluations, matching analytic outputs,
forces, gradient components and per-request cache trace. Full commands and allocation
request are in SUBMISSION.json/run.sbatch. Nested ORCA32CPU timing is diagnostic;
use outer allocation/scheduler time to avoid double-counting total cost.

No new DFT is part of this manifest. More starts, different bounds, alternative
chemistry or a production rescore require a separate recorded scope. Existing
reference mapping preparation and the parent's proposed native candidate checks
are distinct records, not permission to add tasks to this running experiment.
