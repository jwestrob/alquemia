# Explicit three-source execution adapter — experimental, launch held

These commands prepare or inspect real artifacts without molecular calculations.
The source-to-score route is implemented but remains unqualified and unlaunched;
see STATUS.md. No historical canonical member is required.

```bash
THREE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
THREE_PLAN="$PWD/workspaces/pqq_three_source_execution_20260923/preflight_v1/plan.json"
THREE_SCRIPT="$PWD/workspaces/pqq_three_source_execution_20260923/preflight_v1/implementation/pqq_three_source_execution.py"
"$THREE_PY" "$THREE_SCRIPT" dry-run --plan "$THREE_PLAN"
```

Create a new prepared run (new output directory required) using the exact two
existing PLM triples and the existing numerical qualifications:

```bash
"$THREE_PY" scripts/pqq_three_source_execution.py prepare \
  --preparations \
  workspaces/pqq_three_source_20260923/plm_v1/PQQSEQ_07ab500e3df76b30d71c/PREPARATION.json \
  workspaces/pqq_three_source_20260923/plm_v1/PQQSEQ_83440678cbbd658047c9/PREPARATION.json \
  --rank-qualification workspaces/native_gfn2_rank_panel_20260923/run_v2/COMPARISON.json \
  --maxiter-qualification workspaces/adaptive_maxiter_20260922/prepared_v1/collection_1209876.json \
  --agreement diagnostics/pqq_three_source_execution_20260923/PLAN.md \
  --output workspaces/pqq_three_source_execution_example_v1
```

`prepare` only stages coordinates, inputs and immutable implementation. It reuses
exact source preparation, never archived energies or forces. The fixed union
and threefold reference are inherited from each explicitly pinned preparation.
Additional request/preparation details are in
../pqq_three_source_20260923/COMMANDS.md. Actual four-mode selection waits for
fresh paired origin forces. No candidates or calculated values are invented.

Collection/report work on partial or completed output. Use fresh output names;
all six planned sources remain visible when no endpoint has run:

```bash
"$THREE_PY" "$THREE_SCRIPT" collect --plan "$THREE_PLAN" \
  --output workspaces/pqq_three_source_execution_20260923/preflight_v1/inspection_2.json
"$THREE_PY" "$THREE_SCRIPT" report \
  --result workspaces/pqq_three_source_execution_20260923/preflight_v1/inspection_2.json \
  --output workspaces/pqq_three_source_execution_20260923/preflight_v1/inspection_2.md
```

The implemented molecular operation is `execute --plan` with the immutable
script, inside the specified32-CPU/200000-MiB/one-H200 allocation. The wrapper
`run.sbatch` preserves failures and always collects. **Do not submit it at this
checkpoint**: root explicitly held launch after the threefold transfer lost two
reference-group decisions. A later finite agreement/submission receipt must pin
this exact input/method or a separately declared version; no silent threshold or
state change. Existing started runs collect only, without automatic retries.

Focused actual-artifact tests:

```bash
"$THREE_PY" -m unittest discover -s tests -p test_pqq_three_source_execution.py -v
```
