# Masked OMOL response fails the direct-force and curvature gates

All40 real evaluations completed. Numerical derivatives pass; a mechanical
correction remains unavailable. Production baseline and earlier scores are unchanged.

| Check | Passing comparisons | All-case gate |
|---|---:|---|
| Own analytic derivative | 24/24 | True |
| Direct deformation energy versus DFT | 6/24 | False |
| Even/curvature term versus DFT | 19/24 | False |
| DFT gradient plus learned even term | 24/24 | True |

These24 comparisons comprise four representations, two physical directions and
Ca/La/paired-R quantities. They are not24 independent biological observations.
Both GGR representations share one structure; the alpha structures are one
qualified comparison. All are consumed development cases.

## What failed

Maximum direct signed-energy error is0.320106kcal-scale; maximum curvature error
is0.019278. Direct responses fail18/24 comparisons. For example, GGR_extended
La has a learned metal-direction linear change of+0.107640 versus DFT−0.122956
at the declared0.02A step: even the local direction reverses. Five curvature
checks fail: GGR_extended Ca peptide rotation; alpha1F6S metal motion for Ca,
La andR; alpha6IP9 Ca metal motion. No value, direction or threshold was changed.

Anchoring the first derivative to DFT reduces the largest signed prediction
error to0.019487, passing all24 anchored-error checks. That pass uses the
predeclared0.02 absolute floor and does not override the stricter curvature
failures. The earlier POLAR+GB curvature model passed its own screen; the new
masked OMOL curvature has not earned replacement of that model. Neither model
has a validated paired relaxation correction.

The largest own-gradient odd-energy residual is0.000830, within the frozen
numerical criterion. These are actual learned derivatives, not numerical bugs.
The comparison is deliberately against CPCM DFT while OMOL is a masked vacuum
descriptor. That Hamiltonian difference is a limitation; this test does not
isolate which missing physical component causes the disagreement.

## Execution and scope

Job1200888:40successful calls,0failed; eight analytic centers and32signed
energy-only displacements, with byte-identical archived DFT inputs and physical
cap Jacobians. No newDFT, solvent, training, optimization or biological score.
320GPUallocationseconds,5120allocatedcore-seconds,416.572reportedCPU-seconds.
Summed model time20.148532seconds; peakGPU1,031,137,792bytes. Startup dominates
allocation time. Local preparation/report/test receipts are separate.

Two real contract tests passed in7.247s; the scientific integration test was
explicitly skipped before outputs existed and then passed in0.385s after
completion. No fabricated success fixture; all three distinct tests now pass.
Full result: `workspaces/mace_omol_20260917/masked_response_report_v1/result.json`.
Costs: `masked_response_cost_v1.json`; cumulative accounting: V13.

## Next direction

Do not use unanchored masked-MACE relaxation or its failed curvature as a
correction. Continue with the separately declared [static context descriptor](MASKED_SUBTRACTIVE_CONTEXT_PLAN.md): retain DFT local chemistry, subtract the
same learned core from its matched protein. A cached partition prerequisite
can reject that candidate before any new full-protein inference. No production
promotion, broad validation or accuracy gain is claimed.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/masked_response_v1/implementation/mace_omol_response.py report --manifest workspaces/mace_omol_20260917/masked_response_v1/manifest.json --output workspaces/mace_omol_20260917/masked_response_review_report_v1
```

This replays actual receipts; it launches no inference. Use a new output directory.
