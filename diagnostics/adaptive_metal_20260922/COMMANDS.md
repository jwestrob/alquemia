# Prepared joint-metal adapter

No molecular calculation has been launched. Root reviews the physical/numerical
plan and actual continuation before any execution. Existing shared-pool code
does not automatically accept this new protocol; integration requires its own
explicit candidate-union manifest.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
JOINT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
JOINT_RUN="$PWD/workspaces/adaptive_metal_20260922/prepared_v1"
"$JOINT_PY" "$JOINT_RUN/implementation/adaptive_metal_proposals.py" dry-run \
  --manifest "$JOINT_RUN/manifest.json"
```

Preparation already performed (write-once destination):

```bash
"$JOINT_PY" scripts/adaptive_metal_proposals.py prepare \
  --angular-manifest workspaces/adaptive_accommodation_20260922/proposals_v1/manifest.json \
  --agreement diagnostics/adaptive_metal_20260922/PLAN.md \
  --output workspaces/adaptive_metal_20260922/prepared_v1
```

Real-source geometry/parser tests, with no model calls:

```bash
"$JOINT_PY" -m unittest discover -s tests -p test_adaptive_metal_proposals.py -v
```

Read-only collection to a new filename:

```bash
"$JOINT_PY" "$JOINT_RUN/implementation/adaptive_metal_proposals.py" collect \
  --manifest "$JOINT_RUN/manifest.json" \
  --output "$JOINT_RUN/collection_review_v1.json"
```

The executable adapter retains the existing `execute --manifest ...
--pool-collection ...` interface and the existing GPU allocation check, but no
submit wrapper is supplied or run at this preparation checkpoint. All eight
candidate/solvent/score fields in the saved unrun collection remain unavailable.
