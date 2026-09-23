# Opt-in minimal adaptive PQQ command

Run from the repository root, using the pinned CPU environment. Production
`nikasha`/`affordable_workflow` defaults are unchanged. This candidate accepts
explicit canonical PQQ source requests; it does not discover the binding site.

## Fresh source request

The runnable two-crystal example includes actual source paths and hashes, exact
model/chain/metal/PQQ selectors, homologous role selectors and released config:
`diagnostics/pqq_adaptive_candidate_20260923/two_crystal_sources.json`.
For a new compatible source, use the same schema/config and explicitly set its
source pin, unique safe case ID and correct selectors. There is no hidden raw
Ca renaming, assembly selection or unsupported-donor substitution. Predictions
need no biological label. Unknown source-domain evidence remains explicit.

```bash
CANDIDATE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$CANDIDATE_PY" scripts/pqq_adaptive_candidate.py prepare \
  --request diagnostics/pqq_adaptive_candidate_20260923/two_crystal_sources.json \
  --output workspaces/pqq_adaptive_candidate_20260923/my_explicit_run
```

`prepare` and `dry-run` check inputs without molecular calls. `execute` does
fresh source preparation and scoring. Use the generated immutable implementation:

```bash
CANDIDATE_RUN="$PWD/workspaces/pqq_adaptive_candidate_20260923/my_explicit_run"
"$CANDIDATE_PY" "$CANDIDATE_RUN/implementation/pqq_adaptive_candidate.py" dry-run \
  --plan "$CANDIDATE_RUN/plan.json"
sbatch --output="$CANDIDATE_RUN/slurm_%j.out" --error="$CANDIDATE_RUN/slurm_%j.err" \
  diagnostics/pqq_adaptive_candidate_20260923/run.sbatch "$CANDIDATE_RUN/plan.json"
```

The wrapper uses one H200/32CPUs/200000MiB and four concurrent eight-rank ORCA
endpoints. It executes, collects and reports, also collecting when execution
returns nonzero. It submits no further jobs. The already authorized integration
used `two_crystals_v1`; do not rerun it to collect results.

Explicit operations within the allocation / after execution:

```bash
"$CANDIDATE_PY" "$CANDIDATE_RUN/implementation/pqq_adaptive_candidate.py" execute --plan "$CANDIDATE_RUN/plan.json"
"$CANDIDATE_PY" "$CANDIDATE_RUN/implementation/pqq_adaptive_candidate.py" collect --plan "$CANDIDATE_RUN/plan.json" --output "$CANDIDATE_RUN/collected.json"
"$CANDIDATE_PY" "$CANDIDATE_RUN/implementation/pqq_adaptive_candidate.py" report --result "$CANDIDATE_RUN/collected.json" --output "$CANDIDATE_RUN/collected.md"
```

Output names are immutable; choose a new collection name rather than overwriting.
Partial/failed endpoint attempts remain visible and cannot be silently retried.
The exact existing manifest, input hash and scientific state control cache reuse.
A fresh execution requires an allocation; archive replay needs none.

## Exact archive replay

The `archive_replay.json` request pins the two completed recovery pools. It checks
actual source/method/recipe provenance and pool algebra, without calling a model:

```bash
"$CANDIDATE_PY" scripts/pqq_adaptive_candidate.py prepare \
  --request diagnostics/pqq_adaptive_candidate_20260923/archive_replay.json \
  --output workspaces/pqq_adaptive_candidate_20260923/my_archive_replay
"$CANDIDATE_PY" scripts/pqq_adaptive_candidate.py collect \
  --plan workspaces/pqq_adaptive_candidate_20260923/my_archive_replay/plan.json \
  --output workspaces/pqq_adaptive_candidate_20260923/my_archive_replay/result.json
```

Replay is explicitly labeled `exact_archive_replay`; it is not fresh-source
integration. The original released score and its bands stay separate from the
candidate's frozen canonical25 reference. A missing matrix cell produces a
missing candidate result, even when an original score is available.
