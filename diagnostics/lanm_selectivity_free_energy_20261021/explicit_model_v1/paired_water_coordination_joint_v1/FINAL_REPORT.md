# Joint water fit: completed, three checks still fail

The fit stops after22 complete steps with no further measured improvement. Independent
raw-output and objective review passes. Scientific acceptance is **22/25**:9/10
energy and13/15 rigid-force controls. This is a consumed development result.

| Remaining failure | Error | Existing limit |
|---|---:|---:|
| Eight-to-nine-water paired La/Dy energy change | −4.338284 kcal/mol | magnitude≤1 |
| La nine-water rigid-force RMS | 2.042854 kcal/mol/Å | ≤2 |
| Paired expanded-dimer rigid-force RMS | 2.203824 kcal/mol/Å | ≤2 |

The aggregate objective falls1.481382→0.895804. The nine-water energy error improves
only from−4.951230 to−4.338284. Individual nine-water forces improve substantially,
but the paired expanded-dimer force no longer passes its previous limit. An average
improvement therefore does not qualify the model.

All six parameters remain inside their original intervals. At the selected point,
the recorded parameter-response matrix has rank6 and condition149.17. Its local
least-squares projection changes the objective by only1.55e−10; actual final trials
provide no improvement. This supports a local residual limitation of the current
six-parameter response. It does not prove that no other parameter basin can pass,
nor identify one missing physical term.

**Decision:** preserve the failed fit, hold the fresh radial validation and conditional
Cartesian qualification, and investigate additional mechanically meaningful response
parameters before another fit. No automatic continuation, changed thresholds or
combination with the independently running donor fit follows.

The run evaluated338 new parameter vectors through4,120 native programs. Scheduler
allocation:10,776 CPU-seconds for fitting plus9 for the two saved reports. Inputs,
previous failures and the earlier17670 pre-submit-review loss remain preserved.

Evidence: workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/
paired_water_coordination_joint_v1/run_v1/ROOT_FINAL_ADMISSION.json; selected_review_v1/
run_v1/execution/REVIEW.json; residual_review_v1/actual_v1/REPORT.md.
