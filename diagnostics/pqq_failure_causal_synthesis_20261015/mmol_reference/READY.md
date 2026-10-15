# MMOL1770 Ca0: six-cell reference adapter ready

Metadata-only preparation completed; no numerical preflight, molecular call or
submission was performed by this agent. The executable sources are frozen in
`workspaces/pqq_failure_causal_synthesis_20261015/mmol_reference_v1/implementation`.

Manifest SHA256: `10c23bbfcdfacb6a5ff3ea5dbaa6b8661f6ef8cc79681289b4d83b1252b58dc8`.
Inventory SHA256: `18f0fa12759ec0caf8e2a58c561e094d8e920332108e3900855d69469fbb3b35`.
The root [scientific plan](../MMOL_REFERENCE_PLAN.md) is pinned in the manifest.

## Exact staged commands

From the repository root (root owns submission, node placement and actual wakes):

```bash
mmol_work="$PWD/workspaces/pqq_failure_causal_synthesis_20261015/mmol_reference_v1"
mmol_sha=10c23bbfcdfacb6a5ff3ea5dbaa6b8661f6ef8cc79681289b4d83b1252b58dc8
sbatch "$mmol_work/implementation/RUN.sbatch" "$mmol_work/manifest.json" "$mmol_sha" preflight
```

After actual preflight passes, submit the two independent origin scouts:

```bash
sbatch "$mmol_work/implementation/RUN.sbatch" "$mmol_work/manifest.json" "$mmol_sha" execute sample0_Ca_origin_vacuum
sbatch "$mmol_work/implementation/RUN.sbatch" "$mmol_work/manifest.json" "$mmol_sha" execute sample0_La_origin_vacuum
```

Every native action requests 12 task slots on one standard node. One ORCA worker
uses 12 MPI ranks, one numerical thread each, with the original allocation-aware
memory renderer. No memory, walltime, priority or account override. Repeated
scientific work is never inferred from a missing watcher: inspect receipts first.
Per-task locks prevent concurrent writes; the native runner retains valid completed
outputs and refuses incomplete/unverified attempts or stale runtime files.

### Four proposal cells after root's actual scout review

An explicit immutable review JSON must contain:

- `decision: "continue_four_proposals"`;
- `manifest: {path, sha256}` exactly matching this manifest;
- `scout_collections`: the two actual one-cell accepted origin collection pins;
- `interpretation`: root's substantive executable/state/cost review.

Pass its actual path and SHA to each of
`sample0_Ca_adaptive_Ca_vacuum`, `sample0_Ca_adaptive_La_vacuum`,
`sample0_La_adaptive_Ca_vacuum`, `sample0_La_adaptive_La_vacuum`:

```bash
sbatch "$mmol_work/implementation/RUN.sbatch" "$mmol_work/manifest.json" "$mmol_sha" execute TASK_ID --scout-review ABS_REVIEW_JSON --scout-review-sha REVIEW_SHA256
```

This is the predeclared scientific staging gate, not a request for fresh user
permission. Root already owns autonomous continuation within the plan.

### Saved-output comparison

After the desired actual endpoints terminate, pass their actual collection paths
(the script also records specified missing files; it never fills their energies):

```bash
sbatch "$mmol_work/implementation/COMPARE.sbatch" "$mmol_work/manifest.json" "$mmol_sha" ABS_NEW_COMPARISON_OUTPUT_DIR ABS_COLLECTION_1 ABS_COLLECTION_2
```

For the final result supply all six one-cell collections. This job uses one CPU,
performs no molecular evaluation, and writes `RESULT.json` plus `EXECUTION.json`.
Partial common-path results remain available when their required cells exist;
vacuum row minima require all six cells. The fixed candidate tie order is origin,
adaptive_Ca, adaptive_La. Historical composite-selected geometries are reported
separately from vacuum minima. No class band, fitted threshold or source1 endpoint
enters this comparison.

## Narrow implementation review

- Exactly one source, six byte-copied XYZ files and unchanged r2SCAN-3c vacuum
  inputs; both proposal admissions and same-state Ca/La maps must match their
  original pinned task manifest. The inventory-only source1 is never a task.
- Native executor, runtime renderer, scalar/population parsers, and corrected
  printed-TRAH threshold handling are copied unchanged. `COUNTED_COLLECT.py`
  copies the original collection function with the manifest atom count replacing
  its hard-coded 192 population rows; declared output count is also manifest-derived.
- Actual native auxiliary basis, printed singlet and complete population-map
  partition are explicitly checked, alongside the inherited native checks.
  The current protocol does not turn DIIS/SOSCF energy-only endings into gradient
  qualification; solver summaries and state warnings remain in every result.
- Slurm preflight replays the inherited three real output fixtures and all eight
  completed ExaA collections through the count-aware collector, checking exact
  energies, populations, group sums and final solver records against their saved
  results. It computes MMOL all-electron/ECP inventories from actual XYZ there;
  the archived GFN valence count is not reused as a DFT electron count.
- AST parsing, shell syntax, exact snapshot bytes, task membership and existing
  receipt statuses were checked locally. These are metadata/static checks only.
  Scientific preflight is **UNRUN** until root's actual Slurm result is available.
- Curie's independent static source/state and comparison reviews found no blocker. Comparison
  algebra is `E(q)-E(origin)`, then Ca minus La; missing cells never become zeros
  or a smaller candidate pool. Raw native receipt intervals are reported with
  their scope; they exclude queue/preflight/collection and are not total latency.
