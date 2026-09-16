# Global electrostatic challenger: completed, not promoted

**All four quantum endpoints, four ESP checks and 25 distinct TABI tasks
completed. The frozen model fails both partition consistency and numerical
acceptance. Baseline/default remains unchanged; the conditional ten-endpoint
accuracy trial did not run.** [Compact machine-readable result](RESULT.json).

## Main result

Moving Asp303 between the classical and quantum representations changes the
Ca-minus-La descriptor by approximately 18 kcal/mol at every frozen level:

| Numerical setting | qm36 minus qm33, kcal/mol | Frozen limit | Result |
|---|---:|---:|---|
| Primary | +18.075764323 | 2 | Fail |
| Refined surface | +18.039361675 | 2 | Fail |
| Refined tree evaluation | +18.075763683 | 2 | Fail |

The matched archived CPCM boundary-core comparison shifts by +2.204092696
kcal/mol. That is a separate baseline-derived protocol, not canonical fixed-core
PQQ calibration or an inherited challenger threshold. [Archive audit](MATCHED_BASELINE_COMPARISON.md).

At primary settings, the partition change decomposes into vacuum QM
**+61.738603357**, direct core/environment coupling **−43.031559106**, and
whole-protein reaction field **−0.631279928** kcal/mol. The global cavity did
not remove the incomplete cancellation of quantum and direct-field terms.
This does not uniquely identify caps, charge representation or density response.

## Numerical and state checks

| Prescribed check | Observed change/error, kcal/mol | Limit | Result |
|---|---:|---:|---|
| Surface refinement, qm33 paired contrast | −1.396760397 | 0.5 | Fail |
| Surface refinement, qm36 paired contrast | −1.433163045 | 0.5 | Fail |
| Rigid rotation, qm33 paired contrast | +1.713579558 | 0.5 | Fail |
| Tree refinement, largest paired change | 0.000002804 | 0.5 | Pass |
| Translation, absolute paired change | 0.000000093 | 0.5 | Pass |
| Exact repeat | 0 | 0.01 | Pass |
| Native versus independent Coulomb cross accounting, largest error | 0.000000010 | 0.01 | Pass |
| Both isolated algebraic reductions | 0 | 0.01 | Pass |

All five required common-mesh groups have identical actual mesh and source
hashes within their groups. Paired coordinates, protonation, charge ownership,
water inventory, source geometry and full physical cavity checks pass. All
25 solver outputs converged and passed strict native-control/receipt parsing.
These execution successes do not turn the scientific gate into a pass.

Individual reaction-field energies change by 56.34–57.79 kcal/mol on surface
refinement, with partial cancellation in the La/Ca contrast. Rotation changes
individual energies by 5.90/7.62 kcal/mol. The much smaller tree-setting effect
helps locate the numerical problem, but does not uniquely separate surface
discretization from solver convergence. No settings were tuned after inspection.

Full components, endpoint shifts and all 51 checks are in the final assessment.
The two isolated controls check reduction to gas QM plus isolated reaction
energy; they are not a zero-solvation or zero-transfer claim.

## Model and scope

Protocol: `vacuum_r2scan3c_mbis_global_tabi_electrostatic_v1`.

`Etilde_M = E_QM,vac,M + C(q_M,Q) + G_RF,D(q_M+Q)`;
`R_global = Etilde_Ca − Etilde_La`.

Native ORCA 6.1.1 r2SCAN-3c/DefGrid3 gas endpoints, MBIS, fixed ff19SB permanent
protein charges, source-only shared SES cavity; dielectric 1/78.54, zero salt,
298.15 K, 1.4 Å probe, common 1.8 Å metal radius. No CPCM, old isolated-PB
counterterm, extra full Coulomb energy, aquo offset, universal zero or calibrated
classification. This is a frozen-density electrostatic descriptor.

All four new MBIS exterior-potential checks pass (relative RMS 0.01856–0.05323).
That limited quality test does not validate every core/environment interaction.
[Boundary audit](BOUNDARY_AUDIT.md), [frozen numerical protocol](NUMERICS.md),
[approved scope](AGREEMENT.md), [technical recovery](TECHNICAL_RECOVERY.md).

## Execution, failures and measured cost

25 final solver results came from 46 attempts: 25 completed plus 21 preserved
partial attempts interrupted during a persistent low-clock hardware observation.
No quantum retry was needed. Twelve array elements cancelled while pending ran
no calculation. No node settings, queue priority or shared permissions changed.
The original parser rejection of the ORCA SMD credit was fixed by recollection.
[Execution evidence](PERFORMANCE_OBSERVATION.md).

| Allocation | Terminal state | Wall seconds | Allocated CPUs | Allocated CPU-seconds |
|---|---|---:|---:|---:|
| 1199949 | COMPLETED | 932 | 64 | 59,648 |
| 1199952 | COMPLETED | 16 | 4 | 64 |
| 1199956 | COMPLETED | 2643 | 2 | 5,286 |
| 1199959 | CANCELLED by 601 | 1364 | 23 | 31,372 |
| 1199974 | COMPLETED | 4832 | 12 | 57,984 |
| 1199964, executed indices 0–8 | COMPLETED | Per-task receipts | 1 each | 37,584 |

**Total recorded development allocation: 191,938 CPU-seconds
(53.316111 CPU-hours); zero GPU use.** This includes the
interrupted batch and idle portions of multiworker allocations. Sum only unique
parent jobs, never add batch/extern steps again. Slurm counts logical CPUs;
recovery placement verified distinct physical cores.

The primary qm36 pair took 1,713.032 and 1,734.551 seconds on one CPU each
(3,447.583 summed solver wall-seconds). The initial qm33 pair shared a physical
core and took 2,636.471/2,641.909 seconds. Four quantum/MBIS endpoints used
59,648 allocated CPU-seconds. Population analysis itself took 574.688–747.686
seconds per endpoint. These hardware/scheduling differences prohibit a controlled
production overhead ratio. Initial/cached environment preparation and software
validation are not fully included; their absence is explicit, not zero cost.
[Baseline cost audit](COST_REFERENCE_AUDIT.md).

## Tests and judgment

**40 software/real-output integration tests passed, zero skipped, zero failed**
(115.233 s unittest; 115.872317 s process). Scientific calculations are
listed separately above; no fake solver output stands in for a calculation.
[Tests and exact receipts](TESTS.md).

- **Numerically credible at the declared precision? No.** Surface/rotation gates
  fail, although arithmetic, tree refinement and deterministic replay pass.
- **Useful for prediction? Not established.** Physical rejection prevented the
  conditional accuracy trial. No new biological accuracy result follows.
- **Affordable? Partly measured, not demonstrated against a matched baseline.**
  Primary surface solves are modest in CPU count but take roughly half an hour
  each; charge extraction is substantial. The failed model does not earn rollout.

**Recommendation: retain baseline and reject this frozen challenger for
promotion.** This does not reject every global environmental representation.
Jacob subsequently authorized contained new ideas; the separately versioned
[density/field diagnostic](../density_embedding_20260916/REPORT.md) investigates
specific approximation errors without changing this result.

## Reproduce the report

All final artifacts are under `workspaces/global_electrostatic_20260916/`:
`surfaces_v1/{collection_final_v1,assessment_final_v1}.json`,
`surfaces_v1/physical_report_final_v1.md`, `cost_final_v1.json`, and
`software_tests_final_v1/`. Earlier primary/control checkpoints remain intact.
The [runbook](RUNBOOK.md) contains collection commands. A next read-only operation:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/global_electrostatic_assess.py report \
  --assessment workspaces/global_electrostatic_20260916/surfaces_v1/assessment_final_v1.json \
  --output workspaces/global_electrostatic_20260916/surfaces_v1/physical_report_operator_v1.md
```
