# Preparation and actual fixture checks

Run from the repository root. Existing `prepared_v1` is immutable; use an explicit
new output directory for a deliberate repeat. These commands run no energies.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/lanm_global_occupancy_prepare.py \
  --agreement diagnostics/lanm_global_occupancy_20260923/PLAN.md \
  --output workspaces/lanm_global_occupancy_20260923/prepared_v1

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  diagnostics/lanm_global_occupancy_20260923/preparation/test_preparation.py
```

Parent owns molecular feasibility/execution. This preparation module provides no
energy runner, implicit submission or fallback. Source/state files retain exact
parameter/forcefield/source pins. The nine prepared states have 18 physical
endpoints, separate physical multiplicities and native f-in-core singlets.
