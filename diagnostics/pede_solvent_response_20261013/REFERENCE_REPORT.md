# PedE solvent response: independent check complete

The solvent sensitivity persists under independent DFT, but its source-dependent
Ca−La contrast is **+10.03 kcal/mol versus +21.26** with corrected GFN2/ALPB.
The cheap model gives much larger individual Ca and La responses, partly hidden
by cancellation. This supports investigating its electronic/solvent response;
it does not establish that DFT is exact or that the classifier is repaired.

## Interpretation and next decision

The [scientific comparison](comparison/REPORT.md) retains all eight endpoint
energies, both per-metal responses, population changes and the exact algebra.
The [electronic audit](solvent_audit/REFERENCE_ELECTRONIC_LIMITS.md) records
frontier-occupation/gap and stopping warnings. Different electronic and solvent
models changed together, so this test cannot uniquely blame ALPB, the density,
the cavity or a particular fragment. Current state/geometry audits have not
established a shared wrong-only fold or preparation defect.

Proceed with the [fixed-charge native ALPB replay](FIXED_CHARGE_REPLAY_PLAN.md).
It uses archived charge distributions to distinguish frozen-distribution
sensitivity from additional redistribution, including intrinsic compensation.
No new SCC, DFT, search, parameter fitting or classifier promotion is needed.
A failed native-component reproduction blocks interpretation. Production,
original failures and frozen classification bands remain unchanged.

## Actual execution and cost

| Job | Work | Wall seconds | CPUs | Allocated CPU-seconds |
|---|---|---:|---:|---:|
|3354|sample1 La water scout|9580|12|114960|
|3371|sample1 La vacuum|11435|12|137220|
|3372|sample1 Ca water/vacuum|50971|12|611652|
|3373|sample3 La water/vacuum|8052|12|96624|
|3374|sample3 Ca water/vacuum|9529|12|114348|

Total: **8 new molecular endpoints, zero molecular retries,1,074,804 allocated
CPU-seconds, zero GPU-seconds.** Saved audits/reader failures3352/3353/3356 add7;
parser recovery3370 adds2; partial test3376 and final collector3377 add1 each.
The full current phase therefore used1,074,815 allocated CPU-seconds. Allocation
is not a utilization measurement or a future screening-cost estimate.

All molecular workers used12 MPI ranks on node-64-768g-5, a48000MiB actual
memory grant and runtime MaxCore3100MB/rank. Final3372.batch sampled MaxRSS was
17453524KiB; this is Slurm's sampled metric, not a claim of summed process memory.
Its vacuum calculation took12h5m4.199s inside the14h9m31s pair job. The unchanged
native algorithm required843 reported SCF cycles including TRAH work. Printed
MaxIter125 was not a total-work cap; the earlier description is corrected.

The original scout batch failed after successful chemistry because of a parser
lookup. Recovery3370 passed15 parser/real-artifact checks without rerunning it;
3376 passed six partial-comparison fixture checks. Actual molecular integration
comprises the eight separately retained ORCA executions, not those parser tests.
All eight pass the unchanged diagnostic admission policy. No stability search,
new spin/state, additional functional or seeded retry was executed. The prepared
same-state seed contingency is closed because its failure gate was not met.

## Reproduce or resume

Authoritative result:
`workspaces/pede_solvent_response_20261013/comparison_v1/execution/RESULT.json`

SHA256: `0df3271eb35b617bac1c79caf66be24cf5fe52f3da9c631f21dabd7bdf01bba0`.
Exact immutable inputs/executor pins are in [comparison/EVIDENCE.json](comparison/EVIDENCE.json).
Execution receipts and wake acknowledgements are linked from [EXECUTION.json](EXECUTION.json).
[CURRENT.md](CURRENT.md) states the next task; no molecular resubmission is needed.

The already-completed collection command is recorded in comparison/READY.json
and workspace comparison_v1/FINAL_SUBMISSION.json. To inspect the saved result
without recomputation from the repository root:

```bash
cat diagnostics/pede_solvent_response_20261013/comparison/REPORT.md
```

The replay branch will supply its separate runnable Slurm command after interface
review. At this report's creation its build and eight component calls are unrun.
These two repeatedly inspected PedE sources are development diagnostics, not
blind validation or evidence for protein affinity/physiological metal use.
