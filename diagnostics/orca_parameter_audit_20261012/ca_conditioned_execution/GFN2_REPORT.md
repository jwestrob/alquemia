# Ca-conditioned solvent stage: prepared, not executed

The adapter joins the finite CPU MACE batch to the already qualified native
xTB6.7.1/tblite0.7.0 solvent-component path. It reuses the existing worker and
score algebra; it does not recalibrate or introduce new geometry searches.

## Scope and accounting

- All 125 Ca-conditioned sources remain visible: 106 new sources, three actual
  corrected pools reused, 15 original PQQ exclusions and one MMOL2048 exclusion.
- There are 25 declared biological groups and 22 preparation-supported groups.
  C5AXV8, Q88JH5 and Q9Z4J7 each lack all five Ca-conditioned sources.
- The existing corrected canonical bands, origin/adaptive_Ca/adaptive_La order,
  and 0.1 kcal/mol operational origin-retention rule remain fixed. Mathematical
  minima and independently available origin scores are reported separately.
- The existing 1e-12 Å coordinate-copy tolerance applies to paired alias cells;
  residuals are retained and input coordinates are never rewritten by the adapter.
- Exact scientific inputs are deduplicated. At most 1,272 fresh inputs per backend
  are possible before deduplication: two metals × three candidates × two media ×
  106 sources. A task executes native xTB then independent tblite, each one-thread.
  Root must release 14 threads after the MACE batch is terminal before execution.
- Accepted solvent components require the unchanged 1e-6 Eh ALPB-minus-vacuum
  agreement, 1e-5 e per-atom population agreement, convergence and state checks.
  Absolute energy errors remain reported; no estimated offset is subtracted.
- Missing or failed required cells keep accommodation unavailable. A valid origin
  is a separately named result. Interrupted tasks are retained without automatic
  repetition; unstarted tasks may resume. Costs include available failed-call
  receipts and identify started calls lacking terminal receipts.

## Tests on actual artifacts

`REPLAY_TEST.py` normalizes the three completed Ca probe pools into the new matrix
contract, reparses their 36 actual inputs per backend, and reproduces all nine
policy scores and their decisions exactly. Both the corrected summaries and their
raw execution/output pins are checked. These are parser/algebra tests, not new
scientific evaluations or independent biological observations.

The three archived accommodated calls remain two correct and one wrong; origins
remain two inconclusive and one wrong. A deliberately missing-cell copy of one
real fixture confirms that origin survives while accommodation is unavailable.
A changed-charge metadata copy is rejected. The replay workspace cannot execute
scientific workers. No invented energies, structures or labels were used.

The first parser attempt exposed an old MACE receipt schema difference: some
origin receipts store their model pin in the task manifest. The adapter now
validates that actual manifest, matching the existing reevaluation reader.
Historical failed test artifacts remain. The latest final test evidence is
`workspaces/orca_parameter_audit_20261012/ca_conditioned_gfn_parser_v4/RESULT.json`.

The actual 106-source MACE batch and this solvent stage are separate. No solvent
matrix has been prepared or executed before the actual terminal MACE collection.
The final 125-source denominator and fresh-call runtime will be verified then.

## Exact operations

Run from the repository root. These consume the existing finite scientific scope
and the actual terminal MACE collection; they do not restart MACE.

```bash
python diagnostics/orca_parameter_audit_20261012/ca_conditioned_execution/GFN2.py prepare --scope workspaces/orca_parameter_audit_20261012/ca_conditioned_execution_v1/SCIENTIFIC_SCOPE.json --collection workspaces/orca_parameter_audit_20261012/cpu_proposal_batch_preflight_v1/FINAL_COLLECTION.json --output workspaces/orca_parameter_audit_20261012/ca_conditioned_gfn2_v1
python diagnostics/orca_parameter_audit_20261012/ca_conditioned_execution/GFN2.py dry-run --directory workspaces/orca_parameter_audit_20261012/ca_conditioned_gfn2_v1
python diagnostics/orca_parameter_audit_20261012/ca_conditioned_execution/GFN2.py execute --directory workspaces/orca_parameter_audit_20261012/ca_conditioned_gfn2_v1 --workers 14
```

Execution writes per-input receipts, attempt start/finish records and an immutable
`FINAL_COLLECTION.json`. Later recollection requires an unused explicit filename:

```bash
python diagnostics/orca_parameter_audit_20261012/ca_conditioned_execution/GFN2.py collect --directory workspaces/orca_parameter_audit_20261012/ca_conditioned_gfn2_v1 --output-name MANUAL_COLLECTION_001.json
```

Results retain original/accommodated contrasts, selected states, accommodation
works, corrected decisions, exclusions, common-coverage comparisons, numerical
qualification evidence and recorded costs. No new Ca-triple aggregation,
production promotion or biological affinity claim is introduced.
