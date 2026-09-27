# Same-direction environmental response survives this partition change

Four/four real expanded-region endpoints completed. The Thr159 hydroxyl rotation
favors Ca relative to La in both representations: delta=(LaB-LaA)-(CaB-CaA)
is+0.352795817kcal/mol for the original54-atom core and+0.520767257 for the
63-atom core containing complete Thr159 side chain. The change is+0.167971439.
This is a finite representation sensitivity, not a catastrophic discontinuity
or reversal on this case. It does not prove partition convergence or affinity.

| Electronic component | Original core | Expanded core |
|---|---:|---:|
| Ca B-A,kcal/mol | -3.328668903 | -2.588555709 |
| La B-A,kcal/mol | -2.975873086 | -2.067788452 |
| La response minus Ca response,kcal/mol | +0.352795817 | +0.520767257 |
| A positive-arc Ca force,kcal/mol/A | 17.226351367 | 12.210889396 |
| A positive-arc La force,kcal/mol/A | 15.246885300 | 9.575348619 |
| B positive-arc Ca force,kcal/mol/A | 17.772848594 | 15.250980951 |
| B positive-arc La force,kcal/mol/A | 16.049437009 | 12.415776254 |

Individual response/force components change substantially because the expanded
calculation includes Thr internal electronic energy and a new link/boundary
representation. Metal-independent omitted MM terms cancel in the double difference
only if genuinely identical. The differential Ca-minus-La arc force is positive
at both ends in both models; this is consistent with the observed response sign,
not a complete-force or equilibrium claim. Polarization, charge redistribution,
link representation and local quantum chemistry are not separately identified.

Physical source/protons/waters and A/B perturbation are preserved. Expanded
A/B point-charge files are identical; HG1 moves at QM index57. Native original
DefGrid3/TightSCF used for both partitions. The inherited long hydroxyl bond
(1.18537A) remains a limitation. Fine rigid qualification failed on the small
region; the expanded region has no independent numerical error bound. Do not
promote these values to a qualified scorer or treat0.01kcal as a proven bound.
No unknown protein label, biological prediction or optimized geometry is used.

Job1219460:452s,344CPUs,155488allocatedCPU-s,zeroGPU; collector1219461:4s/1CPU.
Five real-fixture harness tests passed. Actual molecular outputs all completed;
no expanded-grid/finite-difference or full-hybrid integration test ran. Original
ML checkpoint/interface and complete hybrid cross terms remain unavailable.
Cumulative molecular workers:42outputs,864128allocatedCPU-s,zeroGPU. This is
one-time development cost, not per-site production throughput.

Next justified direction: check availability of matched preparations for4MAE
and one consumed non-PQQ control under the same frozen hydroxyl selection rule.
That preparation task is independent of numerical qualification; no further
native solver sweep, new classifier or LanM relaxation is authorized by this result.
See RESULT_1219460.json and immutable partition_v1/FINAL_COLLECTION.json.
