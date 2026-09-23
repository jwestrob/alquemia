# Exact six-source experimental integration

This integration uses the two existing PLM triples and the actual envelope
reference. It is not a cohort scorer or a biological accuracy test. Preparation
reuses archived protonation; molecular values are fresh. Default/DFT outputs stay
unchanged. The existing script refuses undeclared groups and automatic retries.

Already submitted: job1211626. Exact commands, immutable plan and wrapper pins are
in `workspaces/pqq_three_source_envelope_execution_20260923/run_v1/SUBMISSION.json`.
Do not resubmit it. The wrapper always attempts collection and retains the
original molecular exit status. It writes `RESULT_1211626.json` and a Markdown
report on termination; these are the authoritative result paths when present.

Read-only preflight of the existing plan:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/pqq_three_source_envelope_execution.py dry-run \
  --plan workspaces/pqq_three_source_envelope_execution_20260923/run_v1/plan.json
```

Explicit new preparation command, with a distinct output path (no execution):

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/pqq_three_source_envelope_execution.py prepare \
  --preparations \
    workspaces/pqq_three_source_envelope_20260923/plm_v1/PQQSEQ_07ab500e3df76b30d71c/PREPARATION.json \
    workspaces/pqq_three_source_envelope_20260923/plm_v1/PQQSEQ_83440678cbbd658047c9/PREPARATION.json \
  --strict-qualification workspaces/strict_native_pool_20260923/run_v1/COMPARISON.json \
  --rank-qualification workspaces/native_gfn2_rank_panel_20260923/run_v2/COMPARISON.json \
  --maxiter-qualification workspaces/adaptive_maxiter_20260922/prepared_v1/collection_1209876.json \
  --agreement diagnostics/pqq_three_source_envelope_execution_20260923/PLAN.md \
  --output workspaces/pqq_three_source_envelope_execution_20260923/reproduction_preparation_v1
```

Report-only recollection, if required after the original wrapper has finished
(use a new result path; never launch a second collector while the wrapper runs):

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/pqq_three_source_envelope_execution.py collect \
  --plan workspaces/pqq_three_source_envelope_execution_20260923/run_v1/plan.json \
  --output workspaces/pqq_three_source_envelope_execution_20260923/run_v1/RECOLLECTION_v1.json
```

The `execute --plan <absolute plan path>` operation is provided by the immutable
script copy and used by the submitted wrapper inside its declared allocation.
No alternative executor or settings are needed. The repeated-start marker blocks
scientific reruns; incomplete results remain collectible and explicit.

## Result schema

- `rows[]`: `protein_id`, `case_id`, unknown `source_evidence`,
  `union_origin_R_model_kcal_mol`, `union_origin_components`,
  `variants.{mathematical,operational}` with accommodated R/decision and same-model
  `delta_R_from_union_origin`; full `pool` selections/endpoint works, `matrix`,
  candidate geometry/boundary/optimizer flags, status/reason. Values are model
  kcal/mol, not affinity probabilities or aquo-referenced DFT scores.
- `groups[]`: exact three `members`, origin median/range,
  `variants.{mathematical,operational}` strict three-source median R/decision,
  member values/range, unknown source evidence and physical anchor.
- Top-level: actual plan/reference/profile pins, six-source/two-group denominators,
  execution receipts, collection pins, experimental status and probe-regression
  limitation. A failed required cell makes that source unavailable; a missing
  source invalidates the corresponding strict group summary.

Raw values must stay in separate candidate fields when joined to existing PLM
exports. Never subtract them from historical DFT or another protocol's raw R.

## Actual cost collection

After job1211626 is terminal, capture job and step accounting (do not use `-X`,
which omits the CPU-time-bearing batch step):

```bash
sacct -n -P -j 1211626 \
  --format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,CPUTimeRAW,TotalCPU,ReqMem,NodeList,MaxRSS \
  > workspaces/pqq_three_source_envelope_execution_20260923/run_v1/ACCOUNTING.txt
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  diagnostics/pqq_three_source_envelope_execution_20260923/summarize_costs.py \
  --plan workspaces/pqq_three_source_envelope_execution_20260923/run_v1/plan.json \
  --accounting workspaces/pqq_three_source_envelope_execution_20260923/run_v1/ACCOUNTING.txt \
  --output workspaces/pqq_three_source_envelope_execution_20260923/run_v1/COSTS.json
```

GPU reservations during scalar phases are counted. Nested source/search/worker/
solver timings are explanatory and must not be summed as additional allocations.
