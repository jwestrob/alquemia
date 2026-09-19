# Archive-only response commands

Run from the repository root. Existing outputs are immutable; the commands
below use new replay paths and refuse to overwrite them. They reuse molecular
outputs and perform no DFT/MACE evaluation or Slurm submission.

```bash
OPENBLAS_NUM_THREADS=1 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/response_probe_archive.py prepare --root . --agreement diagnostics/response_probe_20260919/PLAN.md --output workspaces/response_probe_20260919/replay_features_v1.json
OPENBLAS_NUM_THREADS=1 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/response_probe_archive.py evaluate --features workspaces/response_probe_20260919/replay_features_v1.json --output workspaces/response_probe_20260919/replay_evaluation_v1.json
OPENBLAS_NUM_THREADS=1 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/response_probe_archive.py transfer --features workspaces/response_probe_20260919/features_v2.json --evaluation workspaces/response_probe_20260919/evaluation_v2.json --agreement diagnostics/response_probe_20260919/TRANSFER_PLAN.md --output workspaces/response_probe_20260919/replay_transfer_v1.json
OPENBLAS_NUM_THREADS=1 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_response_probe_archive.py -v
```

The completed run uses corrected `features_v2.json` and `evaluation_v2.json`;
these change direct label metadata without changing numerical features or fitted
weights. `transfer_v1.json` reuses those saved fits. Tests use these real frozen
fixtures; they skip explicitly if missing. Original v1 metadata is superseded.
Native executable inference tests are not part of this archive-only experiment.

The source inventory in the script names the six existing immutable collections
used by this bounded study. Input paths/metadata are recorded in the result.
Do not replace these sources or apply this fitted head to new families without a
new protocol/evaluation declaration. There is no scanner/default integration.
