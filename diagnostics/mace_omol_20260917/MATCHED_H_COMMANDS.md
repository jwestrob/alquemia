# Matched H-normalized hybrid operations

From the repository root. Protocol
`masked_omol_matched_normalized_H_vacuum_hybrid_v1`; the
[plan](MATCHED_H_NORMALIZATION_PLAN.md) fixes all inputs and gates.
Jobs1200950(quantum) and1200951(MACE) were submitted once. Inspect their current
state before executing anything that could duplicate their work. Output paths
must be fresh; scientific outputs and partial failures are never overwritten.

## Validate and collect existing work

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/matched_H_quantum_v1/implementation/mace_omol_matched_h.py dry-run-quantum --manifest workspaces/mace_omol_20260917/matched_H_quantum_v1/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/matched_H_mace_v1/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_20260917/matched_H_mace_v1/manifest.json

/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/matched_H_quantum_v1/implementation/mace_omol_matched_h.py collect-quantum --manifest workspaces/mace_omol_20260917/matched_H_quantum_v1/manifest.json --output workspaces/mace_omol_20260917/matched_H_quantum_collection_replay_v1.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/matched_H_mace_v1/implementation/mace_hybrid.py collect --manifest workspaces/mace_omol_20260917/matched_H_mace_v1/manifest.json --output workspaces/mace_omol_20260917/matched_H_mace_collection_replay_v1.json

/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/matched_H_quantum_v1/implementation/mace_omol_matched_h.py report --quantum workspaces/mace_omol_20260917/matched_H_quantum_v1/manifest.json --mace workspaces/mace_omol_20260917/matched_H_mace_v1/manifest.json --output workspaces/mace_omol_20260917/matched_H_report_v1

/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_omol_matched_h.py -v
```

Missing outputs stay explicit. The actual-result integration test is skipped
until the real report exists; software tests do not stand in for quantum/ML
execution. Original-H results are retained separately for the preparation-effect
comparison; none is substituted for a missing normalized endpoint.

## Fresh reproducible preparation

These commands create inputs and manifests without submitting calculations:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_matched_h.py prepare --prepared workspaces/mace_mechanics_20260916/prepared_v2/manifest.json --audit workspaces/mace_omol_20260917/normalized_H_reuse_audit_v1.json --old-report workspaces/mace_omol_20260917/vacuum_hybrid_report_v1/result.json --agreement diagnostics/mace_omol_20260917/MATCHED_H_NORMALIZATION_PLAN.md --output workspaces/mace_omol_20260917/matched_H_prepared_replay_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_matched_h.py prepare-quantum --preparation workspaces/mace_omol_20260917/matched_H_prepared_replay_v1/preparation.json --output workspaces/mace_omol_20260917/matched_H_quantum_replay_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_matched_h.py prepare-mace --preparation workspaces/mace_omol_20260917/matched_H_prepared_replay_v1/preparation.json --output workspaces/mace_omol_20260917/matched_H_mace_replay_v1
```

The recorded submission argv are in each original manifest directory's
`submission.json`. Quantum execution uses `run_matched_h_quantum.sbatch` and
the existing manifested ORCA executor. MACE uses the existing `run_pilot.sbatch`
with the recorded A5000 resource overrides. Both contain finite task manifests,
reuse locks, preserve failures and collect on exit. No application CPU/time
budget is imposed. Larger cores start first in four16-rank quantum slots;
that ordering changes neither their scientific inputs nor acceptance criteria.
