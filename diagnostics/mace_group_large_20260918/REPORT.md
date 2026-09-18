# Larger grouped model: same predictions, higher cost

All **57 evaluations** completed without failure. Large and medium both pass
**17/22** frozen directional comparisons, with exactly the same five failures.
The declared improvement rule fails. Increasing checkpoint capacity did not
improve discrimination on this consumed development panel.

| Evidence stratum | Medium | Large |
|---|---:|---:|
| PQQ crystal direction | 1/1 | 1/1 |
| Both alpha structures against three GGR structures | 6/6 | 6/6 |
| Khoury domain means against GGR | 6/9 | 6/9 |
| Parvalbumin sites against GGR | 4/6 | 4/6 |

RTX remains below every GGR structure. Parvalbumin CD remains below 2FW0 and
2FVY. All five failing margins become more negative. No successful direction
regresses to a failure, and none of the prior failures becomes a success.
The PQQ margin changes from +57.141809 to +51.626235 model kcal.
These 22 comparisons cover nine biological groups, eight with directional or
supporting labels; they are not 22 independent affinity observations. Khoury
and parvalbumin retain their supporting-evidence limitations. Aequorin remains
the ordered, experimentally unlabelled EF1/EF3/EF4 vector.

## Physical and numerical checks

Protocol `intact_POLAR_large_typed31_group_vacuum_compatibility_v1` changes
the checkpoint only; all 55 medium-panel scientific tasks are identical.
The two additional native controls reproduce actual archived large energies
exactly, forces within 7.084e−13 eV/A and densities within 3.043e−14.
All **69 numerical checks** pass, including native identity, charge closure,
rigid transformation and permutation. All three grouping checks still fail
the retained 2-model-kcal tolerance:

| GGR representation | Connected minus primary, model kcal |
|---|---:|
| 1GLG | −2.8720505231176503 |
| 2FW0 | −2.7798515065805987 |
| 2FVY | −3.055355966615025 |

Qualified comparison count remains zero. No threshold, absolute aquo reference,
PQQ calibration, response correction or production classifier is supplied.
This is a charge-constrained learned vacuum descriptor, not an established
binding free energy or self-consistent electronic Hamiltonian.

## Actual cost and decision

Job 1201533 used **2,475 GPU-allocation seconds (41m15s)**, **39,600 allocated
core-seconds**, and 2,578.268 reported CPU-seconds. Inference totals
2,047.6128651425242 seconds; peak GPU allocation 16,012,100,608 bytes, peak
reserved GPU 22,101,884,928 bytes, and worker RSS 2,033,648 KiB.
Medium transfer used 1,236 GPU-allocation seconds. Large includes two extra
native identity controls, so this approximately twofold ratio is the observed
campaign cost, not an exact matched per-site throughput ratio. All physical/group
preparations and native reference computations are historical reuses; local
preflight/tests/reporting are additional and not fully CPU-profiled.

Numerical reproducibility is credible under the performed checks; grouping
sensitivity remains unqualified. Predictive gain is absent. Execution is
minute-scale and feasible, but the extra capacity did not earn its added cost.
**Close this capacity experiment; do not launch a large canonical expansion.**
Retain baseline/default and the active broader MACE research goal.

Full results and paired margins are in
`workspaces/mace_group_large_20260918/report_v1/`, costs in `cost_v1/`.
RESULT.json pins the compact report; COMMANDS.md gives runnable replay.
Actual result tests are recorded in VALIDATION.md. No DFT, solvent, training,
geometry optimization, label changes or threshold fitting occurred.
