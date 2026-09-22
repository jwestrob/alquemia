# Read-only archived-force replay

Use a fresh output filename. These operations do not call MACE or ORCA.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/adaptive_force_diagnostic.py \
  --manifest workspaces/accommodation_nonlinear_20260920/proposals_v1/manifest.json \
  --manifest workspaces/accommodation_nonlinear_20260920/fold_proposals_v1/manifest.json \
  --agreement diagnostics/adaptive_accommodation_20260922/PLAN.md \
  --output workspaces/adaptive_accommodation_20260922/force_projection_replay_v1.json

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/adaptive_composite_diagnostic.py \
  --source workspaces/accommodation_response_20260920/result_v2.json \
  --agreement diagnostics/adaptive_accommodation_20260922/PLAN.md \
  --output workspaces/adaptive_accommodation_20260922/composite_comparison_replay_v1.json

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  -m unittest discover -s tests -p 'test_adaptive_force_diagnostic.py' -v
```

`project(mapping, q, forces)` returns the exact context geometry, normalized
physical-heavy Jacobian, raw gradients, normalized gradients and mode lengths.
It uses the context Jacobian for force projection, including caps. `preview`
accepts common mode definitions/unit physical vectors and both endpoint loads;
it returns primitive IDs and the projection audit. Neither function optimizes
coordinates, estimates stiffness or modifies the scorer.
