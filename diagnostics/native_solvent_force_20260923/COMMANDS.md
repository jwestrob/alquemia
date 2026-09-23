# Fixed native solvent-force operations

These commands parse existing results; they do not rerun chemistry. The approved
finite experiment is32displaced cells ×3stages, plus eight shared origin gradients
already computed by `native_pool_continuation_20260923`.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
FORCE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FORCE_RUN=$PWD/workspaces/native_solvent_force_20260923/run_v1
"$FORCE_PY" scripts/native_solvent_force.py validate --manifest "$FORCE_RUN/stage0/manifest.json"
"$FORCE_PY" -m unittest discover -s tests -p test_native_solvent_force.py -v
"$FORCE_PY" scripts/native_solvent_force.py compare \
  --design "$FORCE_RUN/design.json" \
  --centers workspaces/native_pool_continuation_20260923/run_v1/stage2/collection.json \
  --output "$FORCE_RUN/comparison_replay.json"
```

All output JSON destinations are new files; immutable originals are not replaced.
The final comparison uses the current analysis script, whose hash is retained in
its result. Each molecular stage executed the pinned implementation in run_v1,
with the same fixed predeclared settings.

For an explicitly authorized fresh reproduction, `prepare` accepts explicit
`--inputs`, `--seeds`, `--agreement`, `--output`; no source editing is required.
The exact executed inputs were:

```bash
"$FORCE_PY" scripts/native_solvent_force.py prepare \
  --inputs diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json \
  --seeds diagnostics/native_xtb_restart_20260922/SEED_AVAILABILITY.json \
  --agreement diagnostics/native_solvent_force_20260923/PLAN.md \
  --output workspaces/native_solvent_force_20260923/preparation_replay
```

`run --design ABSOLUTE_DESIGN_PATH` executes the three declared stages, saving
collections and visible dependent failures after each stage. It requires the
64CPU allocation, eight 8-rank workers and fresh directories. The actual submission
was job 1210338 using `run_v1/run.sbatch`, one64CPU128GiB allocation on
node-128-512g-8gpu-1, no GPU. Its wrapper and exact design are pinned in
`run_v1/SUBMISSION.json`. Do not resubmit a completed experiment or duplicate the
eight q0 cells from the other agent's experiment.
