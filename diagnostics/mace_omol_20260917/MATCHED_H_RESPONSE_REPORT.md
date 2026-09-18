# Matched vacuum hybrid response: ordering improves, qualification remains open

Completed 2026-09-18. **Actual alpha/GGR directions improve from 2/4 to 4/4.**
These are consumed development comparisons involving two biological groups:
two alpha-lactalbumin structures and two representations of the same GGR site.
This is useful evidence for the candidate, not broad or prospective validation.
Production and historical references remain unchanged.

## Actual energy result

| Alpha structure | GGR representation | Static margin | Response margin |
|---|---|---:|---:|
| ALPHA_1F6S | GGR_extended | +4.329632 | +13.202491 |
| ALPHA_1F6S | GGR_connected | +4.211135 | +11.178967 |
| ALPHA_6IP9 | GGR_extended | -0.805099 | +10.796095 |
| ALPHA_6IP9 | GGR_connected | -0.923595 | +8.772571 |

Margins are alpha minus GGR, in the hybrid's kcal-equivalent scale; positive
is the expected direction. All four exceed the frozen 0.02 threshold. No new
absolute calibration, aqueous exchange score, or application of PQQ bands.

Protocol: `matched_normalized_H_vacuum_hybrid_bounded_metal_response_v1`.
The fixed energy is H = native vacuum r2SCAN-3c(core) + masked OMOL(full) −
masked OMOL(core). The matching analytic gradient anchors a quadratic with
whole-model metal curvature inside the unchanged 0.20 Å sphere. Reported
changes above use actual native DFT and matching MACE evaluations at the
predeclared positions, not the quadratic predictions. No entropy/solvent term.
Source atoms, normalized H, core memberships, charges, explicit waters and
scaffold coordinates remain fixed within each paired comparison.

## Physical checks and limits

All eight actual energy decreases and prediction-error checks pass. All
interior/tangent gradient magnitude checks pass. Six of eight endpoints pass
all their checks (27/29 individual checks). The alpha 1F6S Ca and GGR extended
Ca positions fail the radial sign: the model slightly overshoots the actual
constrained minima. Their radial gradient dot displacement values are +0.449876
and +2.331702 kcal-equivalent, respectively. Retain both failures.

The actual response partition change is **1.905028**, passing 2.0; the final
hybrid partition change is **2.023524**, narrowly failing the same frozen gate.
All fully qualified scores therefore remain null. Report v2 adds raw partition
and raw ordering fields even when qualification fails; v1 and its source remain
intact. This reporting improvement changes no coordinates, energy, rule or gate.

All 228 initial calls passed numerical checks: archived center scalar energies
agree exactly; coarse/fine matrices are positive; all axis-gradient checks pass.
Maximum predicted-point refinement energy change is 0.0010763 kcal-equivalent;
maximum displacement refinement is 0.0000406 Å. No eigenvalue clamping.

**Separate judgments:** numerical derivatives and native energy predictions are
credible on these inputs; constrained-minimum and partition qualification are
incomplete. Predictive ordering improves on this consumed pair; transfer is
untested. Affordable production throughput is not established by the development
cost. Recommend pursuing transfer while retaining the production baseline.

## Execution, cost and tests

Jobs 1201385–1201390 all completed. 244 new MACE calls and eight new native
analytic DFT endpoints; zero scientific execution failures/retries. Reused eight
DFT centers, archived scalar centers and two whole GGR gradients with matching
inputs/receipts. No numerical DFT gradients, Hessian, optimization or training.

New development cost: **3,882 GPU-allocation seconds, 101,408 allocated
core-seconds, 40,115.518 reported CPU-seconds**. Native DFT validation took 614
wall seconds on 64 CPUs; native MACE validation took 361 wall seconds on one
A5000/16 CPUs. Initial grids account for the remaining 3,521 GPU seconds.
Preparation: 45.859229 wall / 40.475703 CPU seconds. Primary report: 9.044910
wall / 7.881323 CPU seconds. Native MACE peak allocated GPU memory: 12,023,698,944
bytes. Historical reused calculations cost additional resources; other local
inspection/tests/replay work was not fully timed and is not counted as zero.

Four real native/grid/source/corrupted-input regression tests pass in 16.211s,
none skipped. They replay actual endpoint energies/gradients and assert that
improved raw ordering coexists with failed qualification and null calibrated
scores. Earlier three response-preparation and five matched-H regressions passed.
CPU/GPU frozen preflights passed. Scientific integrations actually ran on ORCA
6.1.1 and the pinned MACE environment.

## Artifacts and next operation

Workspace: `workspaces/mace_omol_hybrid_response_20260918/`.
Primary `native_report_v2/result.json` pins the reporter and actual receipts;
`assessment_v1/result.json` pins all fixed predictions;
`complete_cost_v1/result.json` pins scheduler accounting. Initial/previous
reports remain intact. [Commands](MATCHED_H_RESPONSE_COMMANDS.md) give read-only
report replay; [frozen plan](MATCHED_H_RESPONSE_PLAN.md) defines all original
choices. No jobs remain live for this experiment.

Next scientific priority: test the same model on the already consumed GGR
2FW0/2FVY structural replicas, where direct MACE previously failed. Prepare a
separate declared transfer inventory; do not change radius, bands, chemistry or
qualification criteria to rescue this result. The research goal remains open.
