# Nikasha entrypoint — 22 September 2026

Nikasha is the current product name. `scripts/nikasha` delegates to the existing
manifested workflow using the caller's Python interpreter, arguments, working
directory and environment. The standalone repository has no existing package
console-script configuration, so this adds no installation or dependency.

The [current guide](../../docs/NIKASHA.md), agent-guide heading and README heading
introduce the name. Existing invocations, protocol IDs, calibration, schemas,
cache keys, historical results and repository remote are preserved. No production
scoring implementation changed. The integrated PLM manuscript is the delivery
target; this is not a separate PQQ-only paper plan.

## Actual checks

[CHECKS.json](CHECKS.json) records five pairs of commands through the old and new
entrypoints: general help, standard help, baseline help, ensemble help, and
read-only validation of the completed `standard_1H4I_v1` release plan. All ten
commands exited successfully. Validation output is identical; help changes only
the displayed program name and argparse line wrapping. These are CLI smoke
checks, not new molecular integration runs. Scientific calls: **zero**.

```bash
export NIKASHA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$NIKASHA_PY" scripts/nikasha standard --help
"$NIKASHA_PY" scripts/nikasha standard validate \
  --plan workspaces/pqq_fast_release_20260920/standard_1H4I_v1/plan.json
```

Authorization: Jacob's September 22 Nikasha handoff and delegated parent scope
in [PLAN.md](../nikasha_shared_pool_20260922/PLAN.md). No publication, push, remote
rename or global installation was performed. Unrelated pre-existing README edits
are preserved and excluded from this change's scoped commit.
