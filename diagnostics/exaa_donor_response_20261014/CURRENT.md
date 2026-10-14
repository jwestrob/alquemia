# ExaA donor response: La scout agrees; Ca/control matrix running

**Independent vacuum DFT reproduces the target's large La work:** −26.70294
kcal/mol versus MACE−26.78280. This weakens a MACE-only explanation for that
particular response. Ca–La discrimination, matched-control response and physical
electronic-state uniqueness remain unestablished. Read SCOUT_DECISION.md,
SCOUT_RESULT.json and method_review/PAIR_STATE_REVIEW.md. No classifier gain or
production change is claimed.

## Current execution

The two scouts4109/4110 and join4111 completed. Both molecular endpoints pass
all17 declared native checks. Final TRAH residuals8.804145e−7/5.466705e−7 are
below printed1e−5. Near-final negative-gap/Hessian warnings and level shifting
remain; the final orbital table is not independent stability qualification.

Four of eight cells now complete. The six-cell continuation remains unchanged:

| Job | Task | Last observed status |
|---|---|---|
|4117|sample3 Ca origin|RUNNING|
|4118|sample3 Ca adaptive_La|COMPLETED; energy criterion|
|4119|sample0 La origin|COMPLETED; energy criterion, see below|
|4120|sample0 La adaptive_La|RUNNING|
|4121|sample0 Ca origin|RUNNING|
|4122|sample0 Ca adaptive_La|RUNNING|
|4123|complete eight-cell comparison|PENDING(Dependency)|

Each molecular job uses12MPI ranks on shared standard/node-64-768g-5, normal
priority, no memory/time override. Four active molecular jobs use48 of the current60CPU user
allowance;4122 started after4119 freed its slots. Actual grant/memory renderer
is recorded per execution. All four remaining completion/health watchers plus the
collector wake are live; completed4119/4118 have been acknowledged. Do not duplicate jobs or restart successful scouts.


Control La origin4119 ended normally in3670s on12CPUs, with all16 applicable
checks passing. Unlike the scouts, it ended by the native SOSCF energy criterion,
with actual final orbital gradient4.4651e−5 above printed1e−5; no later TRAH
ending supersedes that residual. This follows the frozen solver-specific admission
rule but is weaker numerical qualification. Preserve the limitation in any control
comparison. See CONTROL_ORIGIN_RESULT.json and method_review/CONTROL_ORIGIN_STATE_REVIEW.md.
No control work or Ca–La differential is available from this single new endpoint.

Target Ca proposal4118 also completed,4736s on12CPUs. All16 applicable checks
pass; it ends on SOSCF energy convergence with orbital gradient1.0425e−5 versus
1e−5, so it also lacks uniform final-gradient qualification. Its matching Ca
origin4117 is still progressing. TARGET_CA_PROPOSAL_RESULT.json retains the
raw residuals and costs. Target-only comparison4125 (1CPU) is submitted to join the
four original sample3 cells after4117 ends, using the unchanged tested comparator.
No molecular evaluations are added; missing control cells remain explicit.
Its completion wake is armed; exact rearm command is target_comparison_v1/SUBMISSION.json.

## Scope and next action

Exact192-atom ExaA two-Ca samples0/3; dry oxidized-PQQ3− hypothesis, charge
Ca−3/La−2, singlets,14caps, origin/adaptive_La unchanged. Common physical path,
not the full operational three-candidate minimum. Native r2SCAN-3c vacuum only;
no GFN2 parameters, added solvent, new folds, optimization or calibrated DFT call.

Workspace: `workspaces/exaa_donor_response_20261014/reference_v1/`. ManifestSHA:
`d37733923c3203f41b73d5a3c1789d396e22f8da69e624df2c2d724d0141bcf1`.
CONTINUATION.json links exact submission/rearm commands and watcher PIDs. Each
job writes `attempts/JOB_execute_TASK/COLLECTION.json`. Full actual joined result
will be `workspaces/exaa_donor_response_20261014/full_comparison_v1/execution/RESULT.json`.
The collector is afterany; failed or missing cells stay explicit.

Read live scheduler and receipts on each event. Inspect actual residuals and
ending solver; do not infer success from scheduler alone or rank estimates from
health monitoring as utilization. No automatic retries. Once Ca work is available,
compare D=W_Ca−W_La on the same path; the full matrix adds the matched-control
D and target-minus-control difference. Use that evidence for the next model
decision, not agreement for La alone. Root owns submissions/integration.

## Actual costs and continuity

Origin6026s×12 + proposal6494s×12 =150240 allocatedCPU-s, zeroGPU. Preflight4108
adds72allocatedCPU-s; join4111 recorded0s integer wall/0.077s TotalCPU. Control
origin4119 adds44040allocatedCPU-s (actual TotalCPU12:04:17; batchMaxRSS11506156K).
Target Ca proposal4118 adds56832allocatedCPU-s. Completed molecular total
251112allocatedCPU-s, zeroGPU; ongoing jobs excluded. These are
targeted development reference costs, not routine scanner costs. Both GBWs and
complete solver outputs are preserved. No old source/calibration is overwritten.

Jacob's discretionary authorization covers this continuation. Vault updated;
substantive email accepted by local relay, inbox delivery unconfirmed. Delayed
4108/4109/4110/4111 notifications do not request reruns. Completed scout callbacks
are acknowledged in their workspace receipts.
