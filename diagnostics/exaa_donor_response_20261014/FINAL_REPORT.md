# ExaA: DFT reproduces the source-dependent donor response

**The full eight-cell comparison supports MACE's local vacuum response on both
sources.** DFT gives a target-minus-control differential of +11.757545 kcal/mol,
versus MACE +11.141984. The large target response is therefore not explained by
a MACE-only local-energy exaggeration in this test. This is useful diagnosis,
not an improved classification or proof that the predicted structure is wrong.

## Complete paired comparison

All entries are kcal/mol. W_M is moved minus origin energy for one metal;
D = W_Ca − W_La. Both metals traverse each host's same saved La-generated path.

| Source / contribution | MACE vacuum | DFT vacuum | DFT minus MACE |
|---|---:|---:|---:|
| Control Ca work | −7.460054 | −7.101245 | +0.358808 |
| Control La work | −12.721550 | −12.479166 | +0.242384 |
| Control D | +5.261497 | +5.377921 | +0.116424 |
| Target Ca work | −10.379323 | −9.567475 | +0.811848 |
| Target La work | −26.782804 | −26.702941 | +0.079862 |
| Target D | +16.403481 | +17.135466 | +0.731985 |
| Target D minus control D | +11.141984 | +11.757545 | +0.615561 |

The two source differentials and their contrast agree within the predeclared
1 kcal/mol interpretive resolution; this is not a measured reference uncertainty.
All eight declared cells are available. Original failures from other experiments
remain intact. The table comes from actual Slurm collector4123, not hand-selected
SCF iterations. FINAL_RESULT.json retains unrounded values and exact provenance.

## What we learned

The independent electronic treatment reproduces both the unusually large response
of sample3 and the smaller response of sample0. Prior added packing changed the
operational composite differential little, and the conditional crystal source
had a much smaller response. These are separate tests with different response
definitions; their energies must not be conflated.

Together they narrow two explanations: the tested extra packing does not dominate
this response, and replacing MACE's local vacuum response with this DFT treatment
would not remove it. Source-dependent strain and shared local preparation/state/
environment assumptions remain live explanations. Neither a crystal nor a fold
is automatically the correct physical state for this conditional calculation.

The specific next test asks whether fixing the metal relative to PQQ forces
protein donors to absorb the uneven geometry. Reciprocal transfer of the original
metal position between the two PQQ frames will isolate that coordinate, holding
all host donor configurations and chemical state fixed. It uses eight MACE
vacuum evaluations, with existing physical gates, before considering further
optimization. See ../exaa_metal_position_20261014/PLAN.md. No new DFT is planned.

## Numerical and scientific limits

Five endpoints end in TRAH/NR below their actual printed gradient threshold;
control Ca origin ends on a passing SOSCF gradient criterion. Control La origin
and target Ca proposal instead meet the frozen native energy-convergence policy,
with final gradients4.4651e−5 and1.0425e−5 above1e−5. Both remain weaker numerical
diagnostics. All applicable original collection checks pass; no threshold was
changed retrospectively. Late negative-gap/Hessian and level-shift warnings,
large summary density residuals and lack of independent electronic-state
stability qualification remain explicit in the result and reviews.

Both Hamiltonians use the same finite dry oxidized-PQQ3− representation. This
neither validates solvent response nor establishes experimental cofactor state.
The common path is distinct from the operational three-candidate row minima;
there are no new class calls, affinities, populations or calibrated reference
bands. These two consumed sources are not independent blind biological tests.
The production discriminator, explicit DFT route and historical results remain.

## Actual work, commands and continuity

Eight native r2SCAN-3c vacuum endpoints completed, each with12MPI ranks:
**522696 allocated CPU-seconds**, zeroGPU. Preflight4108 adds72allocatedCPU-s.
Three collectors consumed0.077,0.087 and0.115 actualCPU-seconds. Latest4122 took
4973s (TotalCPU16:20:44, batchMaxRSS12889428K). This is one-time reference cost,
not routine scanner throughput or measured speedup. All outputs and GBWs remain.

Exact preparation/execution/collection commands and scientific pins are retained
in reference/ and workspace task SUBMISSION.json records. The successful complete
comparison command is in full_comparison_v1/SUBMISSION.json; do not rerun it merely
to recover state. Authoritative result:
`workspaces/exaa_donor_response_20261014/full_comparison_v1/execution/RESULT.json`.
All owned reference jobs are terminal. Pending notifications require acknowledgement,
not molecular retries. Root owns the separately named next experiment.
