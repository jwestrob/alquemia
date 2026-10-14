# ExaA donor response: target differential supported; controls finishing

Independent r2SCAN-3c vacuum reproduces the target's large La-directed response:
**+17.135466 versus MACE+16.403481 kcal/mol**, difference+0.731985. Both metals
lower energy along the same saved physical path; La lowers it more. This weakens
a MACE-only motion artifact explanation. No classifier improvement or biological
metal preference follows. Read TARGET_RESPONSE_REPORT.md and TARGET_RESPONSE_RESULT.json.

## Live state

Six of eight original molecular cells are complete. Only control proposals remain:

| Job | Task | Last observed status |
|---|---|---|
|4109/4110|target La origin/proposal|COMPLETED|
|4117/4118|target Ca origin/proposal|COMPLETED|
|4119/4121|control La/Ca origins|COMPLETED|
|4120|control La adaptive_La|RUNNING|
|4122|control Ca adaptive_La|RUNNING|
|4125|four-cell target comparison|COMPLETED|
|4123|complete eight-cell comparison|PENDING(Dependency)|

Each remaining molecular job uses12MPI ranks on shared standard/node-64-768g-5,
normal priority, no memory/time override;24CPU currently active. Both completion/
health watchers and full-collector wake are armed. No successful cell is rerun.
Exact submissions/rearm commands are linked by CONTINUATION.json.

## Qualification and scope

Target La cells and Ca origin meet their final TRAH residual gates. Target Ca
proposal ends by native SOSCF energy convergence with gradient1.0425e−5 above
1e−5; control La origin similarly ends with4.4651e−5. These follow the frozen
solver-specific policy but lack uniform gradient qualification. Control Ca
origin ends on its native gradient criterion at1.4799e−6. Other large summary
residuals, late orbital-gap warnings and state limitations remain visible.
Convergence does not establish stable unique electronic states.

Exact192-atom ExaA two-Ca samples0/3; dry oxidized-PQQ3− hypothesis, charges
Ca−3/La−2, singlets,14caps, origin/adaptive_La unchanged. Common physical path,
not the complete operational candidate pool. Native r2SCAN-3c vacuum only;
no GFN2 parameters, new solvent, folds, optimization or recalibration.
Production and historical source/score records remain unchanged.

## Next precise action

Read actual4120/4122 collections and
`workspaces/exaa_donor_response_20261014/full_comparison_v1/execution/RESULT.json`
when collector4123 completes. Compare individual works, D=W_Ca−W_La, then
D_target−D_control, with solver qualifications retained. Failed/missing cells
stay unavailable. Use the complete contrast to decide the next physical-model
test; no automatic rerun or extra optimization is prepared.

Workspace: `workspaces/exaa_donor_response_20261014/reference_v1/`.
ManifestSHA `d37733923c3203f41b73d5a3c1789d396e22f8da69e624df2c2d724d0141bcf1`.
Task results live in `attempts/JOB_execute_TASK/COLLECTION.json`; raw outputs,
execution receipts and GBWs remain under `tasks/`. Completed notifications have
ROOT_ACK or ROOT_REVIEW receipts. Duplicate events do not authorize reruns.

## Cost and continuity

Six completed molecular cells total379092allocatedCPU-s, zeroGPU; ongoing work
excluded. Preflight4108 adds72CPU-s; joins4111/4125 each report0integer wallseconds
with0.077/0.087actualCPU-seconds. All exact costs and hash pins remain in prior
endpoint/target reports. These are targeted reference costs, not scanner costs.
Jacob's discretionary authorization continues. Repository and vault record the
result; substantive email accepted by local relay (inbox delivery unconfirmed).
Full matrix interpretation remains pending the two control proposals.
