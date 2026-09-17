# Coupled response: curvature transfers, paired scores remain unsupported

All four representation rows pass the frozen local DFT, grid-refinement, analytic short-gradient and normal/tight-SCF checks. These are three structures from two biological groups: GGR and alpha-lactalbumin. Both alpha structures remain one group, with condition-qualified evidence distinct from the GGR same-assay direction. All are consumed development.

| Representation | La quadratic model | Ca quadratic model |
|---|---|---|
| ALPHA_1F6S | trust_region_exceeded | eligible_for_independent_DFT_minimum_validation |
| ALPHA_6IP9 | trust_region_exceeded | trust_region_exceeded |
| GGR_connected | unstable_or_unresolved_curvature | eligible_for_independent_DFT_minimum_validation |
| GGR_extended | unstable_or_unresolved_curvature | eligible_for_independent_DFT_minimum_validation |

## What was learned

The independently checked cheap local curvature transfers from GGR to both alpha structures, including the coupled metal/peptide direction. Maximum even-energy error is0.009271kcal/mol; each comparison passes its predeclared maximum of0.005kcal/mol or25% of the DFT even term. This is not a uniform0.005kcal/mol accuracy claim. Maximum fine/coarse quadratic change on the complete trust square is0.00008718kcal/mol, below the frozen0.01 gate. All gradients retain their original units, physical mappings and signs.

GGR La has negative combined curvature in both representations. Alpha La has positive two-coordinate curvature but predicted optima outside the fixed0.02Angstrom/1degree rectangle. Alpha6IP9 Ca is also outside. These facts reject a relaxation correction under this model; they do not prove a unique biological failure mechanism. Signs of path curvature away from a stationary geometry do not establish whole-protein basin stability or an affinity label. No eigenvalue was clipped and no trust region was enlarged.

Three Ca predictions qualify for independent validation: alpha1F6S and both GGR representations. Their predicted changes are−0.012543,−0.012360 and−0.032826kcal/mol, respectively. The fixed geometries are pinned in minimum_v1 before new outputs; jobs1200771/1200772 completed3DFT+6short calls. No La contribution is available. Even successful Ca validation cannot supply a paired score.

## Accounting and limitations

The model is the frozen PLAN.md expression: native DFT energy/gradient anchors, MACE+GB core curvature, and the whole-minus-core short component with its linear gradient. Exterior atoms are constrained, not integrated out. The physical protein uses original source H to match DFT; known stretched-H geometry remains a limitation. Long-range scaffold electrostatics is omitted. No duplicate full solvation energy or vertical environmental score is added.

Initial grid jobs1200743/1200749/1200750/1200767 all completed:36newDFT,116newMACE,116newGB,106newshort calls;20exact GGR points reused in each core method. No failed scientific attempt. One preparation filename collision stopped before compute; its unmeasured preparation cost is explicit. Successful preparation took57.513s wall/55.865CPU-s. Initial-grid allocated cost:2,465GPU-s and169,040core-s; reported actual CPU122,419.295s. DFT wall2,025s on64CPU. Conditional validation cost is additional; see the completed result below. See costs.json for raw accounting and memory.

Numerical credibility: passed for the tested local approximations and grids. Predictive improvement: none demonstrated; paired response is unavailable for all four rows. Affordability: the small learned components run quickly, but this finite-grid development workload is not an accepted production score. The baseline remains unchanged.

Exact matrices, gradients, source/receipt hashes and excluded minima are retained in the source assessment linked by grid_result.json. COMMANDS.md provides preparation, collection and assessment commands. No calibrated threshold or entropy correction is enabled.


## Completed independent energy-change validation

All three eligible Ca endpoints pass. The predicted positions were fixed before
these outputs; no new geometry selection, gradient or follow-up optimization:

| Endpoint | Predicted change | Actual DFT+J change | Prediction minus actual |
|---|---:|---:|---:|
| ALPHA_1F6S_Ca | -0.012543191 | -0.012645539 | 0.000102349 |
| GGR_connected_Ca | -0.032826295 | -0.032138699 | -0.000687596 |
| GGR_extended_Ca | -0.012359982 | -0.011423404 | -0.000936578 |

Values are kcal/mol. Every actual change lowers the model energy, and all errors
are below0.001kcal/mol. No DFT gradient was calculated at these predicted points;
this validates their energy change, not exact stationarity under DFT+J. The
component terms are retained in minimum_result.json. No La endpoint qualifies,
so all four paired corrections and the paired GGR partition test remain
unavailable. No zero substitution, new threshold, entropy or accuracy claim.

Conditional jobs cost95GPU-s,28,208allocatedcore-s and11,094.124reportedactual
CPU-s. Slurm allocated64CPU to the quantum job despite the48-task request.
Combined initial-grid and conditional-validation cost is2,560GPU-s and197,248
allocatedcore-s, with133,513.419reportedactualCPU-s. These are development costs:
39newDFT,116newMACE-core,112newshort-component and116newGB calls. Earlier engine
validation and cached GGR work have separate receipts and are not erased.

Four real-fixture preparation/result tests pass in0.614s, including exact source/cap
geometry and rejection of a corrupted real prediction that tries to promote an
out-of-trust endpoint. The retained result is useful but narrow: MACE curvature
and the DFT-anchored response construction work on these supported small Ca
motions. This model does not yet improve the paired discriminator. Retain the
baseline; continue the MACE goal using a separately declared candidate.
