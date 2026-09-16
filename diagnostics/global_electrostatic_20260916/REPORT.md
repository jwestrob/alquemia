# Global electrostatic challenger: progress report

**2026-09-16 — primary partition gate failed; remaining numerical/component
checks running.** Moving Asp303's representation changes the Ca-minus-La
contrast by **+18.075764323 kcal/mol**, exceeding the frozen 2 kcal/mol limit.
The baseline/default remains unchanged. The conditional accuracy stage will
not launch under this gate; no predictive-accuracy result is claimed.

## First partition result

Both representations' primary La/Ca pairs converged and collected with valid
receipts. The qm36-minus-qm33 change decomposes as follows:

| Contribution to partition change | kcal/mol |
|---|---:|
| Vacuum quantum contrast | +61.738603357 |
| Direct core/environment Coulomb contrast | -43.031559106 |
| Whole-protein reaction-field contrast | -0.631279928 |
| Total global contrast | **+18.075764323** |
| Archived matched CPCM baseline contrast | +2.204092696 |

The large change is in the incomplete cancellation of quantum and direct-field
terms. This component accounting does not uniquely identify a cap, charge
boundary, density-response or other physical cause. Refinement/rigid/component
checks are still needed to establish numerical behavior. Their settings and
tolerances remain unchanged; all running checks will be retained.

Primary checkpoint: `workspaces/global_electrostatic_20260916/surfaces_v1/assessment_primary_v1.json`
(six computed solver tasks, zero failed collections, nineteen unavailable at
collection). The matched baseline is a separate protocol, not an inherited
threshold or calibration for the challenger.

## What this tests

Protocol `vacuum_r2scan3c_mbis_global_tabi_electrostatic_v1` combines native
vacuum r2SCAN-3c core energy, endpoint MBIS charges interacting with fixed
protein charges, and one whole-protein TABI reaction-field energy. It uses the
same frozen nuclear states and a common source-atom cavity. It is a static
electrostatic descriptor, not a complete binding free energy.

The approved first gate is 1H4I qm33/qm36 × La/Ca: four quantum endpoints and
the [frozen 25-task solver schedule](NUMERICS.md). The conditional five-pair
accuracy trial has **not run**. [Approval](AGREEMENT.md).

## Completed evidence

All four endpoints converged with the intended gas-phase Hamiltonian and
native ECP convention. All four new MBIS electrostatic-potential checks passed.

| Endpoint | QM energy, Hartree | MBIS total charge, e | ESP relative RMS |
|---|---:|---:|---:|
| qm33 La | −1754.987616371759 | −0.999981 | 0.053225 |
| qm33 Ca | −2400.952491961507 | −1.999989 | 0.028505 |
| qm36 La | −1983.471978518628 | −1.999979 | 0.036323 |
| qm36 Ca | −2629.338467381298 | −2.999984 | 0.018564 |

Unrounded values and endpoint receipt/output hashes are in
`workspaces/global_electrostatic_20260916/partition_tasks_v1/collection_verified_v1.json`;
individual ESP receipts are under `esp_v1/<endpoint>/quality.json` in that
workspace. The original collection falsely matched the SMD author-credit text;
parser correction and recollection resolved it without any QM rerun.

The [boundary audit](BOUNDARY_AUDIT.md) verifies common physical coordinates,
charge ownership and the declared cap approximation; its nine real-artifact
integrity/algebra tests passed. These tests are not solver validation.

## Surface checks in progress

Job **1199956** completed the initial two primary solves. Job **1199959** completed
both isolated controls, but its 21 full-protein solves ran on a socket reporting
about 0.9 GHz under full utilization. The healthy initial pair was preserved;
only the affected owned batch was stopped for same-input technical recovery.
Array **1199964**, indices0–8, restarts nine exact tasks, one solver per allocation.
Its twelve remaining elements were cancelled while still pending, with zero
execution, and grouped into the existing twelve-worker runner as **1199974**.
Both request no SMT sharing and exclude that node. Initial observations on the
replacement nodes show distinct physical cores at about 2.8–3.2 GHz. Existing
queue limits remain unchanged. All interrupted attempts and cost are retained.
The master scientific manifest remains the same 25 tasks.

[Measured execution evidence](PERFORMANCE_OBSERVATION.md) distinguishes the
observed low frequency from its unresolved cause; no overheating was detected.
These runtime observations must not be attributed to the electrostatic model.

The first two surfaces built in about four seconds each and have identical
actual mesh hashes. All four primary endpoints and both isolated-core controls
completed; **the other nineteen results remain pending at this checkpoint**.
Meshing or isolated-control success alone does not pass that gate.

The combined [software checks](TESTS.md) ran 37 tests: 36 passed and one was
explicitly skipped pending all 25 solver outputs. These tests launch no
scientific executable. The [technical recovery record](TECHNICAL_RECOVERY.md)
distinguishes parser/provenance repairs from scientific reruns.

The [matched baseline archive](MATCHED_BASELINE_COMPARISON.md) has a
qm36-minus-qm33 contrast difference of +2.204092696 kcal/mol. These are
baseline-derived boundary-test cores, distinct from canonical PQQ calibration.

## Measured cost so far

| Completed work | Job | Wall s | Allocated CPUs | Allocated core-s | Batch peak RSS KiB |
|---|---:|---:|---:|---:|---:|
| Four QM + MBIS endpoints | 1199949 | 932 | 64 | 59,648 | 8,010,252 |
| Four ESP checks | 1199952 | 16 | 4 | 64 | Not captured |
| Interrupted surface batch, including two completed isolated controls | 1199959 | 1,364 | 23 | 31,372 | 2,017,012 |
| Initial primary whole-protein pair | 1199956 | 2,643 | 2 | 5,286 | 158,968 |

Accounting records: `workspaces/global_electrostatic_20260916/{quantum,esp}_accounting.json`.
ESP accounting reports zero RSS; that is not treated as zero memory use.
Preparation and unfinished solver costs are not included in these totals.
Slurm CPU-time units count allocated logical CPUs on the observed SMT nodes;
they are not physical-core counts. The two initial solvers share one physical
core. The recovery requests one solver per core and records actual placement.

The [same-core cost audit](COST_REFERENCE_AUDIT.md) recovers real no-MBIS and
old CPCM+MBIS receipts. Population analysis takes 574.688–747.686 s per current
endpoint, a substantial part of the quantum wall time. Different hardware,
MPI sizes and schedules prevent a controlled overhead or production-cost ratio.

## Three separate judgments

| Question | Current answer |
|---|---|
| Is the model numerically credible? | **Numerical checks pending; primary partition consistency failed.** Quantum/ESP, input and primary solver checks passed; the 18.08 kcal/mol partition shift exceeds the frozen tolerance. |
| Does it improve La/Ca accuracy? | **Untested.** Stage 2 has not run; this partition experiment supplies no new biological validation. |
| Is it affordable? | **Pending.** Quantum/ESP costs are measured; complete solvent/preparation cost and a controlled baseline comparison are unavailable. |

The [accuracy-input audit](ACCURACY_INPUTS.md) recovered all five frozen core
pairs. Exact canonical 1H4I/4MAE environments currently fail terminal templates;
GGR/alpha require peptide/water mapping and explicit environment-state handling.
No missing atom, alternative source or guessed heterogen charge was substituted.

## Final gate and next action

**Gate: failed at primary level; Stage 2 remains unexecuted.** Collect all 25
solver tasks from the master manifest, preserve failures, then assess the
predeclared 0.5 kcal/mol numerical/rigid, 2 kcal/mol partition and 0.01 kcal/mol
accounting tolerances. Review full-environment preparation before any Stage 2
quantum execution. Partial successes cannot trigger it.

Exact collection, assessment and reporting commands are in [RUNBOOK.md](RUNBOOK.md).
The parent executor will replace this checkpoint with the completed gate,
component results, total measured cost and resulting recommendation.
