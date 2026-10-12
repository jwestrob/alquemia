# CPU proposal batch command

Prepared, reviewed and launched by root as worker1620438/watcher1620439.
Do not start a duplicate; the command below documents the active execution and
explicit recovery after it is terminal. The
`preflight_v1` directory is deliberately the execution directory too, so the
106-source physical checks were not repeated just to change its name.

From the repository root, run this one foreground command:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -u workspaces/orca_parameter_audit_20261012/cpu_proposal_batch_preflight_v1/CPU_BATCH.py execute --manifest workspaces/orca_parameter_audit_20261012/cpu_proposal_batch_preflight_v1/manifest.json
```

The exact same command performs an explicit restart: terminal sources and actual
completed origin/cross receipts are reused. A begun incomplete optimizer remains
unavailable; the command does not begin it again. The run lock and recorded
process start identities reject overlapping owned runs. All-terminal restart
rebuilds the collection without loading a model.

The command owns two optimizer parents and two persistent seven-thread CPU
models. Each parent computes only while its own child waits; the parent waits
during inference. Maximum active numerical threads:14. Each new attempt must
pass both Ca/La archived force/energy controls on each worker before any new
source starts. The original optimizer bounds, maxiter200 and stopping tolerance
remain unchanged; there is no additional time or call ceiling.

## Receipts

Workspace: `workspaces/orca_parameter_audit_20261012/cpu_proposal_batch_preflight_v1/`

- `manifest.json`, `DRY_RUN.json`, `CPU_BATCH.py`, `PLAN.md`: frozen preparation.
- `attempts/NNN/lane{0,1}/`: process ownership, qualification, model output,
  parent logs, cost and failure records.
- `attempts/NNN/QUALIFICATION_RELEASE.json`: all four real controls passed.
- `sources/<case>/`: origin receipts, real paired-force selector, original
  optimizer results/evaluations, common-pool cross receipts, terminal source.
- `attempts/NNN/FINAL_COLLECTION.json`: immutable attempt collection.
- `FINAL_COLLECTION.json`: latest collection index for root's existing watcher.

The collection records `terminal:true` only after all106 source results exist;
individual failures can still make scientific pools unavailable. It preserves
`terminal:false` and pending IDs after interruption. Available native pools have
all six admitted MACE cells. Scalar solvent/composite/classification fields stay
unavailable until the separate corrected-GFN2 stage.

No completion watcher is embedded here. Root uses the existing stage watcher
with this exact workspace and `FINAL_COLLECTION.json` result name.
