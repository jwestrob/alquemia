# Coupled donor response does not improve discrimination

**Completed:4/12 raw directions,0/12 qualified comparisons.** All eight
2FW0/2FVY comparisons remain wrong. Every margin is worse than the preceding
metal-only response. This result closes this fixed-path candidate as an accuracy
improvement. The production baseline and earlier records are unchanged.

## Actual matched comparisons

Model kcal-equivalent; each entry is alpha minus GGR. Expected direction is
greater than0.02. All are consumed development, representing two biological
groups. Structural/core replicas are not independent biological observations.

| Alpha | GGR | Source geometry | Metal only | Coupled donors |
|---|---|---:|---:|---:|
| 1F6S | 1GLG_extended | 4.3296 | 13.2025 | 10.7588 |
| 1F6S | 1GLG_connected | 4.2111 | 11.1790 | 8.3254 |
| 1F6S | 2FW0_extended | -16.8310 | -5.7198 | -9.0562 |
| 1F6S | 2FW0_connected | -20.0623 | -11.6526 | -15.4689 |
| 1F6S | 2FVY_extended | -15.1231 | -4.0334 | -8.6146 |
| 1F6S | 2FVY_connected | -16.9500 | -8.9037 | -14.8351 |
| 6IP9 | 1GLG_extended | -0.8051 | 10.7961 | 7.2333 |
| 6IP9 | 1GLG_connected | -0.9236 | 8.7726 | 4.7999 |
| 6IP9 | 2FW0_extended | -21.9657 | -8.1262 | -12.5817 |
| 6IP9 | 2FW0_connected | -25.1971 | -14.0590 | -18.9944 |
| 6IP9 | 2FVY_extended | -20.2578 | -6.4398 | -12.1401 |
| 6IP9 | 2FVY_connected | -22.0848 | -11.3101 | -18.3606 |

The scalar predictions already gave4/12; actual native DFT retains4/12.
The energy changes agree within the frozen tolerances on all16 endpoints,
so replacing these predicted changes by native values does not repair the
failed directions. This does not validate the learned forces or general curvature.

## Physical checks

All16 native energies decrease along the selected motions. Energy prediction
and full-MACE grid replay checks all pass. Ten of16 endpoints pass every native
check;58of64 individual checks pass. The six failures concern all-mode gradient
prediction: both1F6S endpoints,6IP9Ca, and the extended Ca core of each GGR
structure. All connected GGR endpoints pass native checks yet the two additional
structures remain misordered. No gate was relaxed and no gradient-based
uncertainty, entropy, stationary minimum or free energy is inferred.

| GGR | Response difference | Final score difference | Frozen criterion |
|---|---:|---:|---|
| GGR_1GLG | 2.314831 | 2.433327 | Both absolute differences <=2; failed |
| GGR_2FVY | 4.393556 | 6.220537 | Both absolute differences <=2; failed |
| GGR_2FW0 | 3.181417 | 6.412753 | Both absolute differences <=2; failed |

These measure the scoring procedure's sensitivity to the two core
representations, including their different selected paths. They are not a
decomposition of a fixed-coordinate boundary error, and no unique physical
cause is assigned. All qualified scores remain unavailable.

## What ran and cost

Protocol `matched_vacuum_hybrid_coupled_donor_path_v1`. Source-defined metal
translation, donor chi rotations and all coordinating peptide crankshafts;
maximum physical-heavy displacement0.20A. A fixed learned nonlinear full
potential plus the center Cartesian DFT-minus-learned core tangent selects
among four prescribed nonzero points and the reused center. Native evaluation
then computes DFTcore+Tfull-Tcore at that selected point. Same normalized H,
vacuum native r2SCAN-3c, core inventory, waters, assembly and source microstates.
The [frozen plan](COUPLED_PATH_PLAN.md) defines every criterion.

All16 new analytic DFT and96 MACE calls completed; zero scientific failures.
Jobs1201398–1201402,1201407,1201408. Quantum1489wall seconds/64CPU; native
MACE820wall seconds/16CPU/oneA5000. Total **2112GPU-allocation seconds,
129088allocated core-seconds,90201.055reported CPU-seconds**. Learned forward
time1281.373228s; peak GPU allocation12023698432bytes. These totals exclude
the separately recorded archived centers and additional local work.

Local preparation67.267659wall/66.249144CPU-s. Four actual-fixture tests passed
61.767s, zero skips: analytic mixed-coordinate maps, real/corrupted manifests,
actual selection/corrupted point, actual native algebra and complete report
replay. Earlier geometry/projection tests also passed and remain separately
recorded. No test substituted fabricated scientific outputs. Positions-only
normalization avoids unused Jacobians and reproduces the existing geometry
exactly; neither backend nor scientific paths changed.

## Judgment and outputs

Numerical energy replay is credible; force/representation qualification is
incomplete and predictive usefulness did not improve. Production cost is not
established on matched hardware, and this candidate has not earned additional
routine compute. **Retain the baseline; stop expanding this path candidate.**
The physical coordinate and gradient interface remains useful research support.

[Compact result](COUPLED_PATH_RESULT.json); [operations](COUPLED_PATH_COMMANDS.md).
Full source pins, actual endpoints, components, receipts and cost record are in
`workspaces/mace_site_response_20260918/coupled_report_v1` and
`coupled_complete_cost_v1`. [Comparison PDF](../../workspaces/mace_site_response_20260918/coupled_figure_v2/comparison.pdf),
[SVG](../../workspaces/mace_site_response_20260918/coupled_figure_v2/comparison.svg).

Next read-only replay:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_site_path_native.py report --preparation workspaces/mace_site_response_20260918/coupled_native_v1/preparation.json --output workspaces/mace_site_response_20260918/coupled_report_replay_v1
```
