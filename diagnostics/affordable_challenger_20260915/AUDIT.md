# Live audit and delivery status — 2026-09-15

**Recommendation: retain baseline.** Implementation and preparation checks are
available; the environmental scientific pilot is queued, not validated.

## Current and historical state

The checkout started on `main` at
`bc79e0c84c2db8fff7619ace5b36a47af3be1d62`, with substantial existing edits and
untracked scientific implementations. No reset, stash, baseline edit, queue
reprioritization, production rescore or push was performed. Only the new
affordable-development files and a scoped coordination entry belong to this
delivery. [Implementation inventory](IMPLEMENTATION_REVIEW.json) records the
live dependencies and their preserved workspace copies. Public HEAD alone
does not reproduce those previously uncommitted implementations.

| Area | Finding / disposition |
|---|---|
| Current routing and scorers | Left untouched; generic/PQQ typed31 v2 and their reference policies remain accessible. |
| PQQ fixed-core v3 | Preserved calibrated reference experiment; all 27 released input/XYZ/output/receipt sets verify, and parsing/algebra reproduce published values. |
| Calibration | 25/25 controls, 8.602932 kcal/mol gap. Motif, core charge and composition already separate this panel; incremental DFT value is not established by its accuracy alone. |
| Crystal transfer | 1H4I 7.517843566 and 4MAE 38.089733042 kcal/mol on the released display gauge. These cases are consumed development data now. |
| Prior point-charge embedding | Approximately 64–68 kcal/mol partition discontinuity rejected that tested implementation. The result is retained, not generalized to all environment models. |
| Generic peptide donor | Confirmed C/O + two H defect. New version retains the actual amide N and appropriate caps; canonical fixed-core PQQ did not contain this defect. |
| Non-PQQ panel | Six already prepared sites repaired without energetic rescoring: GGR, ordered aequorin EF1/EF3/EF4, carp parvalbumin CD/EF. Original preparations remain intact. |
| External controls | At audit time the six non-PQQ and 1KB0 endpoint results were absent. Their energies were not opened to develop this challenger. No calbindin opposing-site labels were promoted. |
| Legacy work | BVS, 12-6-4, FEP, restrained minimization and QM/MM are existing work, not new proposals. Historical within-lanthanide failures do not by themselves test La/Ca discrimination. The audited HansLanM QM/MM status was stale; no missing final quantum result was inferred from a basin ensemble or prose. |
| APBS backend | Official 3.4.1 Linux release installed only in the task workspace; executable/version checked. No scientific APBS solve has run yet. |
| Environmental preparation | Four qm33/qm36 endpoint states reproduce exactly; common physical cavity and endpoint-invariant protein charges verified. Both fixed-core endpoint environments fail on incomplete terminal Lys, with no added atom/cap or guessed charge. |
| Response research | Analytic extraction and cap chain-rule mapping implemented; coordinate Jacobians tested on a real bond. No real gradient fixture or validated curvature backend; numerical mechanical/entropy corrections remain unavailable. |

The earlier `environment_preparation_failures.json` is an immutable attempt
record. Its four qm33/qm36 identifier errors were fixed by normalizing blank
insertion-code identifiers. [Verified preparation](environment_preparation_verified.json)
supersedes those four statuses; the two terminal-Lys failures remain current.

## Protocols and score policy

| Protocol ID | Status |
|---|---|
| `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3` | Existing reference, unchanged |
| `generic_peptide_amide_vertical_native_r2scan3c_v3` | New preparation repair, six sites; no calibration or new energies |
| `native_r2scan3c_cpcm_mbis_endpoints_v1` | New opt-in property-extraction pilot; same native quantum method and coordinates |
| `frozen_mbis_apbs_transfer_descriptor_v1` | New opt-in transfer descriptor; scientific integration pending, uncalibrated |

Released bands come directly from the calibration result, not a refit:
Ca-supported raw R ≤ −405404.9238187139; La-supported raw R ≥
−405396.32088671124 kcal/mol. In its optional aquo reporting gauge these are
14.857129202922806 and 23.460061205609236. Values between the bands are
indeterminate. The conversion factor is exactly the released 627.509474
kcal/mol/Hartree. The gauge is not a new reference-compatibility claim.
No new protocol inherits these bands. See [energy accounting and assumptions](MODEL.md).

## Tests and actual cost

[Verification receipt](verification.json): **19 tests, 17 passed, 2 scientific
integration tests explicitly skipped** while jobs are pending. Tests use real
archived/prepared artifacts; malformed-input tests identify their corrupted
copies. They cover published extraction/signs/bands, missing references,
nonconvergence, paired coordinates/charges, proline, overlap, numbering and
insertion codes, incompatible altlocs, missing/broken peptide neighbors,
physical cap Jacobians, cache settings, finite pilot and budget rejection.

No new ORCA endpoint, APBS energy, analytic gradient or displacement DFT has
completed. Environmental identity, convergence, transform and partition
checks therefore remain **unrun**, not passed. Full failure/retry scientific
integration remains untested; partial ORCA output is retained and requires an
explicit fresh retry within the shared budget. Low-level automatic replay is
disabled; retained partial results cannot be overwritten or counted as success.

Measured reference cost: job 1198050 used **154 s × 344 allocated CPUs =
52,976 core-seconds for four concurrent endpoints**, with batch MaxRSS
12,401,728 KiB. This is not a matched per-site production ratio. The measured
repeat environment preparation took **19.56792014092207 wall seconds**, with
19.913344 user CPU seconds, 0.348653 system CPU seconds and process peak RSS
207,776 KiB (CPU/RSS include imports). Earlier preparation/dependency work was
not fully instrumented and is not assigned zero cost. No new GPU allocation.

Scientific judgments are separate: **numerical credibility pending;
incremental predictive information untested; ordinary-score affordability
unestablished**. A canonical-panel success or a favorable partition sign
would not establish broad validation.

## Running / next action

- ORCA job **1198934**: six endpoints, maximum eight including retries;
  high-level admission budget 158,928 allocated core-seconds.
- Dependent solver job **1198939**: at most six ESP checks and 24 six-solve
  APBS states (currently four supported states imply 18 variants), budget
  52,976 allocated core-seconds. Total development ceiling 211,904.
- Both pending at this audit; Slurm estimated the first start at
  **2026-09-17 00:51:50 local cluster time**, subject to change.
- Runtime budgets stop additional admission; local policy prohibits
  subprocess timeouts. In-flight overruns remain recorded and disqualify cost
  acceptance. Missing interrupted-job cost blocks further admission.
- Task-owned watchers record terminal `sacct` receipts. The solver writes
  `workspaces/affordable_challenger_20260915/solver/REPORT.md` on completion;
  a startup/infrastructure failure remains visible in its Slurm stderr.

Next command: `squeue -j 1198934,1198939`. See [runnable operations](OPERATIONS.md)
and the [current paired ledger](comparison_reviewed/REPORT.md).

The work remains open at scientific integration: inspect actual endpoint/ESP
outputs, collect physical checks and accounting, then assess the frozen
challenger. Do not consume new controls, tune parameters after seeing results,
rerun successful endpoints, or promote the model automatically.
