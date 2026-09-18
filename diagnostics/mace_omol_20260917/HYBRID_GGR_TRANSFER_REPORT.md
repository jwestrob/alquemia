# Matched hybrid response fails the additional GGR structures

Completed 2026-09-18. **None of the eight new alpha/GGR orderings is corrected.**
Across all three GGR structures and both core representations, actual response
improves the directional count from2/12 to4/12, exclusively through the previously
inspected1GLG cases. The all-structure robustness requirement fails. All cases
are consumed development and represent two biological groups, not12independent
affinity observations. Production and historical references remain unchanged.

## All actual comparisons

Positive alpha-minus-GGR is the declared expected direction. Units are the
hybrid's kcal-equivalent scale; these are actual native energies at the fixed
predicted positions, not quadratic estimates or aqueous binding free energies.

| GGR structure/core | 1F6S static | 1F6S response | 6IP9 static | 6IP9 response |
|---|---:|---:|---:|---:|
| GGR_1GLG_extended | +4.329632 | +13.202491 | -0.805099 | +10.796095 |
| GGR_1GLG_connected | +4.211135 | +11.178967 | -0.923595 | +8.772571 |
| GGR_2FW0_extended | -16.831014 | -5.719836 | -21.965744 | -8.126231 |
| GGR_2FW0_connected | -20.062349 | -11.652577 | -25.197080 | -14.058973 |
| GGR_2FVY_extended | -15.123068 | -4.033402 | -20.257799 | -6.439798 |
| GGR_2FVY_connected | -16.950049 | -8.903720 | -22.084780 | -11.310115 |

The eight new margins improve by8.05–13.84kcal-equivalent, but remain negative
(−4.03 to−14.06). No averaging, favorable structure/core selection, label change,
threshold fitting or displacement-radius increase was used. Four1GLG comparisons
are reused exactly from MATCHED_H_RESPONSE_REPORT.md. The same separate
GGR direct Ca-favoring and alpha condition-qualified La-favoring evidence applies.
Sugar/crystallization states and pH7 versus assay qualifications remain explicit.

## Physical checks

All156initial MACE calls, eight initial DFT endpoints,16native MACE gradients
and eight native DFT endpoints completed, with no scientific execution failure.
All four whole center energies exactly reproduce scalar archives. All24odd-axis
gradient checks pass; fine/coarse curvature is positive and refinement passes
(maximum predicted-point energy change0.0021538kcal-equivalent). All eight
native positions retain donor membership and all nonmetal coordinates.

Six of eight native endpoint models pass all checks (21/25 individual checks).
Both extended-core Ca states fail energy prediction. 2FVYextendedCa also fails
the interior gradient criterion;2FW0extendedCa fails the boundary radial sign.
Their prediction errors are−1.473190 and−1.342042kcal-equivalent. All actual
energy changes are downward; a decrease alone does not qualify the model.
The four connected-core endpoints pass energy and interior-gradient checks,
yet their alpha/GGR directions still fail. Fixing a numerical or gradient
qualification issue alone therefore cannot establish predictive success.

| Structure | Response partition change | Final hybrid partition change | Frozen gate |
|---|---:|---:|---|
| GGR_2FVY | 3.043337 | 4.870317 | Both fail <=2.0 |
| GGR_2FW0 | 2.701406 | 5.932741 | Both fail <=2.0 |

All qualified scores remain null. Original1GLG failures remain in the reference
report; no successful transfer is substituted for an absent qualification.

## What this teaches us

Metal displacement supplies a real, directionally helpful correction but does
not remove source-structure sensitivity or the dependence on the QM boundary.
Changing only the metal position is insufficient on these inputs. The connected
and extended representations disagree about the size of the Ca response even
though they share the same whole-system curvature. This identifies a limitation
of this assembled correction; it does not uniquely identify solvent, caps,
protonation, scaffold compliance or a specific physical cause.

**Numerical credibility:** scalar replay and derivative/grid calculations pass;
the cheap response approximation fails the two extended Ca checks and partition
qualification. **Scientific usefulness:** margins improve, but no new direction
is corrected and broader robustness fails. **Affordability:** development cost
is measured below; a production cost advantage has not been established.

**Recommendation:** abandon this metal-only correction as the route to a broad
replacement, retaining its gradient/mapping machinery and the production
baseline. Do not spend another pilot merely tuning its radius/threshold or
optimizing its Hessian implementation. The goal remains open; the next model
choice should address a different limitation and be declared before new outputs.
No further scientific pilot is launched in this report.

## Implementation and measured cost

New source-backed transfer and conditional-native modules reuse the established
OMOL worker, native ORCA recipe and runners. Same protocol:
`matched_normalized_H_vacuum_hybrid_bounded_metal_response_v1`;
new transfer stage `matched_hybrid_GGR_transfer`. Preparation reconstructs both
existing58/111atom GGR policies and reproduces archived normalized1GLGinputs.
No production code path, scorer default, core rule or old output is migrated.

New transfer compute: **3,758 GPU-allocation seconds; 172,768 allocated core-seconds; 107,344.786 reported CPU-seconds**.
Jobs1201391–1201394,1201396–1201397 completed. The two64CPUquantum batches took
887and873wall seconds; GPU centers/grids/native took267/1513/1513/465seconds.
16newnativeDFT and172MACEcalls, zero scientific retries; historical reference
calculations are additional. Source-preparationV1's filename-only failure is
preserved;V2 corrected that technical path before inference. Local unprofiled
preparation failure, manifest generation, tests and housekeeping are not zero.
Full precise timings and memory are in complete_cost_v1/result.json.

Two real source/mapping tests passed7.613s. The finite-manifest/corrupted-charge
test initially skipped before manifests existed, then passed16.204s and16.158s
with the final fixed-descriptor guards. Actual native-position/mapping/corrupted-
prediction test passed5.568s. Final actual-output report regression replays all
12contrasts, source energies and unavailable qualified scores in5.216s. No fake energies or successful scientific outputs were used.
Frozen native CPU/GPU preflights passed; unavailable integrations were not
silently replaced. No own jobs remain.

## Reproduction and records

[Plan](HYBRID_GGR_TRANSFER_PLAN.md), [exact commands](HYBRID_GGR_TRANSFER_COMMANDS.md).
Workspace `workspaces/mace_omol_hybrid_transfer_20260918/`: config.json,
prepared_v2,initial_v1,assessment_v1,minimum_v1,report_v1,complete_cost_v1.
The primary result pins its reporter, all actual receipts and the original
reference report. Interim static_report_v1 remains intact. The next runnable
operation is read-only report replay using the documented report command and
fresh `report_replay_v1` output. No new inference is needed to review this result.
