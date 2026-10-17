# Module 2 saved comparison

Root owns binding and Slurm execution. `COMPARE.py` has no molecular backend,
optimizer, polling or submission path. Its only non-standard dependency is NumPy,
used to load the actual pinned LBT gradient files during Slurm analysis.

## Bind actual results

```bash
python diagnostics/lanm_coupled_mechanics_20261017/comparison/COMPARE.py \
  --stage template --repo "$PWD" --output /ABS/OWNED/INPUTS.json
```

Fill six ordered `searches[].result` values with exact paths or pins. Keep
unavailable cells as `null`; list earlier failed result paths in
`attempt_results`. Bind the original successful `scout_result` once at top level,
the final `lbt_pool_result`/`mex_pool_result`, and snapshots of actual completion
receipts in `cost_receipts`. An omitted receipt remains explicitly missing.
Do not include watcher-only jobs or another session's jobs in these cost records.

Actual completed four LBT results reported by root are:

- `workspaces/lanm_coupled_mechanics_20261017/lbt_run_7CCO_La_v1/execution/RESULT.json`
- `workspaces/lanm_coupled_mechanics_20261017/lbt_run_7CCO_Dy_v2/execution/RESULT.json`
- `workspaces/lanm_coupled_mechanics_20261017/lbt_run_7CCN_La_v2/execution/RESULT.json`
- `workspaces/lanm_coupled_mechanics_20261017/lbt_run_7CCN_Dy_v2/execution/RESULT.json`

The original scout is `workspaces/lanm_coupled_mechanics_20261017/lbt_scout_v1/execution/RESULT.json`.
The three failed LBT v1 attempts and two failed Mex v1 attempts remain prior
attempts. Mex v2 and the paired pools must be bound from actual completed files,
never presumed from submission. Read CURRENT and actual receipts first.

```bash
python diagnostics/lanm_coupled_mechanics_20261017/comparison/COMPARE.py \
  --stage bind --config /ABS/OWNED/INPUTS.json --output /ABS/OWNED/CONFIG.json
```

Template/bind only hash/read metadata and can run on the orchestration host.
Run the following **inside root's normal one-CPU Slurm allocation**, using the
existing MACE environment Python so NumPy is available:

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python \
  diagnostics/lanm_coupled_mechanics_20261017/comparison/COMPARE.py \
  --stage compare --config /ABS/OWNED/CONFIG.json --output /ABS/OWNED/execution
```

The output directory must not exist. `RESULT.json` is preserved even on analysis
failure. Complete report generation does not mean every scientific cell is
available or that any search converged. Native work, exact lineage, force norms,
missing-menu suppression and accounting checks are part of this single saved-data
analysis; no additional numerical test suite or chemistry is needed.
