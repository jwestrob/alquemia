# Finite-dielectric PCM recovery remains unsuccessful

All seven allocated groups completed their declared inventory. Fresh native
Model objects removed the error-state contamination between endpoints, but
the scientific/numerical qualification still fails: **312 of 350 checks pass**.
All 18 new forward roles started; 12 failed to converge within 1,200 iterations.
Six new roles succeeded, including the two zero-source controls. Two earlier
successful coarse states were reused exactly. Missing energies remain null.

The available complete 2FW0 coarse pair has R=-43.0844134162 kcal/mol;
the refined pair has R=-38.9020125954. No complete 2FVY pair is available at
any resolution. No between-structure comparison, affinity score or class can
be reported for this model. The coarse reciprocity failure also remains.

This concludes the declared finite-PCM recovery. Increasing the iteration
allowance did not yield a reliable inexpensive backend. The separate conductor
model is a different physical approximation and has its own qualification.
No wider finite-PCM computation is proposed from this result.

## Actual executions and cost

Jobs 1201286–1201292 used 17,340 summed allocation-wall seconds and
1,109,760 allocated core-seconds. Slurm reports 1,073,219 CPU-seconds at
whole-second precision. These are summed independent jobs, not elapsed time
on a single node. No GPU, DFT, MACE, force or biological scoring calls.
All failures are included. Reused earlier runs retain their separate costs.
Read-only final reporting took 32.496134 seconds and made no solver calls.

Manifest SHA256:
`bac125e5580cf3af1c84c102876007be20d27d501a4a0e0e97059a1ce74c24aa`.
Artifacts under `workspaces/mace_omol_20260917/`:

- `ddx_source_recovery_v1/collection.json`: all actual role receipts.
- `ddx_source_recovery_report_v1/result.json`: frozen numerical checks.
- `ddx_recovery_cost_v1/{sacct.txt,result.json}`: allocation/CPU accounting.

Baseline/default, geometry, charges, labels and tolerances remain unchanged.
