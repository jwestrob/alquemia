# Matched hybrid GGR transfer operations

Run from the repository root; use a fresh output path for each replay. The
frozen plan is HYBRID_GGR_TRANSFER_PLAN.md. Inspect jobs before retrying.

## Preparation

Actual configuration and source pins:
`workspaces/mace_omol_hybrid_transfer_20260918/config.json`.
The original prepared_v1 filename-only failure is preserved. Actual complete
preparation is prepared_v2, and initial manifests are in initial_v1.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_hybrid_transfer.py prepare --config workspaces/mace_omol_hybrid_transfer_20260918/config.json --output workspaces/mace_omol_hybrid_transfer_20260918/prepared_replay_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_hybrid_transfer.py initial --preparation workspaces/mace_omol_hybrid_transfer_20260918/prepared_v2/preparation.json --output workspaces/mace_omol_hybrid_transfer_20260918/initial_replay_v1
```

## Existing runners and actual executions

Jobs1201391centers,1201392grid2FVY,1201393grid2FW0,1201394quantum.
Each manifest directory's submission.json stores the actual executable argv,
manifest hash, job receipt and submission time. The established GPU runner uses
one A5000/16CPU/64474MiB; quantum uses64CPU in four16-rank slots. No duplicate
executor, queue reprioritization or changed allocation policy.

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_hybrid_transfer_20260918/initial_v1/centers/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_hybrid_transfer_20260918/initial_v1/centers/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_hybrid_transfer_20260918/initial_v1/quantum/implementation/mace_omol_matched_h.py dry-run-quantum --manifest workspaces/mace_omol_hybrid_transfer_20260918/initial_v1/quantum/manifest.json
```

Runner writes collection_job_JOB.json on exit. Current read-only collect and
collect_quantum commands take --manifest. They retain missing/nonconverged
states and accepted actual receipts. Initial assessment requires all156MACE
and8nativeDFT endpoints and the original numerical gates:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_hybrid_transfer.py assess --initial workspaces/mace_omol_hybrid_transfer_20260918/initial_v1/initial.json --output workspaces/mace_omol_hybrid_transfer_20260918/assessment_v1
```

## Conditional fixed native points

Once the assessment exists, prepare its eligible points without changing the
radius, donor inventory or criteria:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_hybrid_transfer_minimum.py prepare --assessment workspaces/mace_omol_hybrid_transfer_20260918/assessment_v1/result.json --output workspaces/mace_omol_hybrid_transfer_20260918/minimum_v1
```

This creates at most8nativeDFT and16MACE gradients; no trajectory or optimization.
The same two existing runners consume its quantum/mace manifests. After their
completion, the report preserves all12raw comparisons and qualification limits:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_hybrid_transfer_minimum.py report --prepared workspaces/mace_omol_hybrid_transfer_20260918/minimum_v1/preparation.json --quantum workspaces/mace_omol_hybrid_transfer_20260918/minimum_v1/quantum/manifest.json --mace workspaces/mace_omol_hybrid_transfer_20260918/minimum_v1/mace/manifest.json --output workspaces/mace_omol_hybrid_transfer_20260918/report_v1
```

As of this checkpoint only initial jobs have run. Native preparation/report
commands are implemented but await actual eligible predictions; do not describe
those integrations as completed. Production and old experiments remain unchanged.
