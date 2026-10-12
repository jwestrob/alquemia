# Ca-conditioned transfer: finite scope frozen

The next test covers every saved Ca-conditioned reference fold under the same
expanded own-three policy. The scientific scope is prepared; this document
launches no molecular work.

- 125 declared sources from 25 biological groups.
- 106 new admitted sources; three completed corrected pools reused exactly.
- 15 prior PQQ exclusions and one MMOL2048 preparation exclusion retained.
- 22 biological groups have supported sources; the complete Ca source sets for
  C5AXV8, Q88JH5 and Q9Z4J7 remain unavailable.

`SCIENTIFIC_SCOPE.json` pins all source/state/map identities, current model and
optimizer/pool settings, actual reused rows, and corrected canonical bands.
No calibration is repeated. The prior A0AC Ca sample4 error remains in the
reused denominator. No new Ca-triple aggregation rule is introduced.

Final scope:
`workspaces/orca_parameter_audit_20261012/ca_conditioned_execution_v1/SCIENTIFIC_SCOPE.json`

SHA256: `6e31585165092f10a8294d31c7def13331fe0e6196949b5a3eb94dc75cf59ccb`

Heisenberg owns the actual CPU runner, two persistent seven-thread workers,
qualification barrier and restart handling. Root owns launch and the subsequent
qualified GFN2 component stage. The agreed MACE terminal collection is
`workspaces/orca_parameter_audit_20261012/cpu_proposal_batch_v1/FINAL_COLLECTION.json`.
No downstream solvent implementation is a prerequisite for launching the
reviewed MACE batch.

## Cost forecast

The actual CPU scout used 33 force calls on one 202-atom source. If all 106 new
sources required the same 28 optimizer calls plus four origin/cross calls,
plus four initial worker qualification calls, that would be 3,396 forwards.
The measured eight-thread per-call rate implies about **5.83 hours of ideal
two-lane inference time**, plus optimizer and orchestration overhead. This is
an extrapolation, not a promised finish time or compute cap: the planned workers
use seven threads, new sources range from169–215 atoms, and optimizer work varies.
The downstream maximum is 1,272 GFN2 inputs per backend before exact deduplication.
Actual calls, failures, CPU/wall time and RSS must come from execution receipts.

## Reproducible preparation

The already completed scope dry-run verified every new origin/map/state pin,
all three reused qualified pools, full denominator and frozen reference identity.
The preparation operation uses only standard-library artifact reads:

```bash
python diagnostics/orca_parameter_audit_20261012/ca_conditioned_execution/PREPARE_SCOPE.py --config diagnostics/orca_parameter_audit_20261012/ca_conditioned_execution/CONFIG.json --output workspaces/orca_parameter_audit_20261012/ca_conditioned_execution_replay_v1
```

Use the existingv1 scope for the authorized batch. A new output path is shown
only for reproducibility; no completed chemistry needs repeating.
