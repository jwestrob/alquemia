# Opt-in ten-fold request, preparation and archive replay

This is a preparation/replay interface. It does not submit or execute new
molecular calculations. Existing production commands and defaults are unchanged.
The supported group is exactly five Ca-conditioned plus five La-conditioned
sources, including its actual archived canonical La member; selected crystals are
singletons. Three-fold membership is a separately named future experiment.

## Inspect the completed real checks

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
UNION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$UNION_PY" scripts/pqq_union_candidate.py dry-run \
  --plan workspaces/pqq_union_candidate_20260923/two_groups_v3/plan.json
"$UNION_PY" scripts/pqq_union_candidate.py dry-run \
  --plan workspaces/pqq_union_candidate_20260923/crystals_v2/plan.json
"$UNION_PY" -m unittest discover -s tests -p test_pqq_union_candidate.py -v
```

Actual per-source tables and strict aggregates:

- `workspaces/pqq_union_candidate_20260923/two_groups_v3/RESULT.json` and `REPORT.md`.
- `workspaces/pqq_union_candidate_20260923/crystals_v2/RESULT.json` and `REPORT.md`.

## Build an explicit request and replay into a new directory

The archive bundle pins actual source preparations, the unchanged reference,
released static comparison and completed union/adaptive matrices. All input paths
and group IDs are parameters; no source-code editing is needed. The example below
uses existing A0A3 and Q9 data and zero new molecular calls.

```bash
"$UNION_PY" scripts/pqq_union_candidate.py request \
  --sources diagnostics/accommodation_controls_20260920/PQQ_ALL250_SOURCES.json \
  --groups a0a3f2yly8-pqq-la_model q9z4j7-pqq-la_model \
  --archive workspaces/pqq_union_candidate_20260923/ARCHIVE_v1.json \
  --output workspaces/pqq_union_candidate_20260923/review_request.json
"$UNION_PY" scripts/pqq_union_candidate.py prepare \
  --request workspaces/pqq_union_candidate_20260923/review_request.json \
  --source-mode reuse-exact \
  --output workspaces/pqq_union_candidate_20260923/review_v1
"$UNION_PY" scripts/pqq_union_candidate.py replay \
  --plan workspaces/pqq_union_candidate_20260923/review_v1/plan.json \
  --output workspaces/pqq_union_candidate_20260923/review_v1/RESULT.json
"$UNION_PY" scripts/pqq_union_candidate.py report \
  --result workspaces/pqq_union_candidate_20260923/review_v1/RESULT.json \
  --output workspaces/pqq_union_candidate_20260923/review_v1/REPORT.md
```

Output directories/records are immutable. Use a new explicit output path for a
new preparation. `--source-mode fresh` uses the existing normalization,
Ca-residue-only adapter and protonation policy instead of cached source
preparation; it still performs no energy/force calculation and retains failures.
Fresh scientific scoring is unavailable in this interface pending a qualified
execution profile. The dry-run's future counts are declarations, not executions.
An unknown source never acquires a same-named archived score.

## Existing crystal selectors and rebuilding the archive bundle

```bash
"$UNION_PY" scripts/pqq_union_candidate.py request \
  --sources diagnostics/accommodation_controls_20260920/PQQ_ALL250_SOURCES.json \
  --crystal-sources diagnostics/pqq_fast_release_20260920/SOURCES.json \
  --crystals 1H4I 4MAE 1KB0 \
  --archive workspaces/pqq_union_candidate_20260923/ARCHIVE_v1.json \
  --output workspaces/pqq_union_candidate_20260923/review_crystal_request.json
"$UNION_PY" scripts/pqq_union_candidate.py archive \
  --preparation workspaces/consistent_context_20260922/prepared_v1/PREPARATION.json \
  --crystals workspaces/consistent_context_20260922/crystal_controls_v1/CRYSTALS.json \
  --comparison workspaces/union_adaptive_20260923/transfer225_v1/COMPARISON_v1.json \
  --reference workspaces/union_adaptive_20260923/CANONICAL_REFERENCE_v1.json \
  --output workspaces/pqq_union_candidate_20260923/review_archive.json
```

Preparation first fixes a fragment union from every supported declared source,
then reconstructs each source's own coordinates/caps and compares physical state
identities across members. It does not force different folds to have equal
coordinates. All excluded members remain in the declaration and strict summaries.
The old numerical profile stays explicit; cheaper stopping/rank profiles and
three-source union calibration are not silently substituted.

Canonical membership is checked against the pinned historical source record; it
is not a user-assigned state anchor. This interface is limited to existing
reference groups and explicit crystal controls. It does not advertise arbitrary
unknown ten-fold PLM proteins as qualified inputs. A general-input adapter must
separate state anchoring from calibration membership and La4 exclusions.
