# 7CCN-La completion attempt: further model lowering, still nonconverged

Run15458 completed the declared **80 additional accepted steps and 80 distinct
evaluations**, reaching **172 accepted steps / 172 evaluations** cumulatively.
It terminated at `accepted_step_limit`, with `converged=false`.

For the exact step92→172 segment, the saved model composite decreases by
**15.585663 kcal/mol**. The fixed-FF ALPB contribution decreases while the
whole-MACE vacuum contribution increases. **This segment has no independent
reference validation.** The previous PBE0 comparison covered step80→92 only.

## Saved producer results

| Incremental work, step172 minus step92 | kcal/mol |
|---|---:|
| Whole-MACE vacuum | +5.308810 |
| Fixed-FF ALPB | −20.894474 |
| Composite | −15.585663 |

These are the producer's saved component differences. Within this model,
the solvent contribution outweighs the vacuum-energy increase. This does not
establish the accuracy of either contribution or of the new endpoint.

| Full256 composite force summary | Step92 | Step172 |
|---|---:|---:|
| Maximum atomic force, kcal/mol/Å | 5.674648 | 4.824511 |
| RMS atomic-vector force, kcal/mol/Å | 1.086248 | 1.045572 |

The final full-force L2 norm is **16.729158 kcal/mol/Å**. Both convergence
requirements remain unmet: maximum atomic force ≤0.1 and RMS atomic-vector force
≤0.03 kcal/mol/Å. These endpoint summaries do not establish monotonic force
improvement, force accuracy or a minimum.

The original source was not rebased. All saved trace and trial geometry gates
pass, including the original bond, chirality, peptide-basin, severe-overlap and
state checks. All 256 physical atoms remained movable. Passing these guards
does not establish electronic-state or force qualification; both qualification
flags remain false.

## Scope and limitations

This was one La search on the consumed, Lu-conditioned 7CCN peptide. It made
**80 MACE calls and 80 ALPB calls**, with no DFT, GFN, MM, opposite-metal or other
search calls. It reused the prior checkpoint and retained the existing source
and chemistry assumptions.

The prior step80→92 reference supported that earlier incremental energy
lowering under its stated SCF/density and state caveats. It cannot validate
step92→172. Conditional terminal chemistry, the vacuum model with post-added
fixed-FF solvent, nonconverged geometry and unresolved physical-state accuracy
remain limitations. No metal preference, affinity, general accuracy or minimum
claim follows. The producer records `no_automatic_followup=true`; root owns
any subsequent scientific decision.

## Actual cost and evidence

| Allocation | Elapsed seconds | CPUs | GPUs | Allocated CPU-seconds | Allocated GPU-seconds |
|---|---:|---:|---:|---:|---:|
| Preflight15457 | 13 | 1 | 0 | 13 | 0 |
| Run15458 | 58 | 12 | 1 | 696 | 58 |
| **This attempt total** | | | | **709** | **58** |

These are scheduler allocation costs, not utilization. Prior searches and
reference allocations are excluded. The producer reports 56.389137 seconds
application elapsed time, including 6.878603 seconds model loading,
10.254591 seconds in MACE calls and 4.221440 seconds in ALPB calls.
[Exact accounting](run_result_v1/COST.json) preserves the raw scheduler records.

The [lossless archive and SHA manifest](run_result_v1/EVIDENCE_PINS.json) retains
all **80 evaluations**, full gradients and native component receipts,
**82 checkpoint files**, the full **81-row trace**, saved parent admissions,
source/configuration, launch records and logs. No original artifact was replaced.

- [RESULT.json](run_result_v1/execution/RESULT.json.gz), uncompressed SHA256:
  `9159ea43cf76692fb9ad40acbf199c1a6b146acdf8dc21ff68a045350a764a10`.
- [SEARCH.json](run_result_v1/execution/SEARCH.json.gz), uncompressed SHA256:
  `4d59031becf63f38381a33733bbda62fef4ebe48e8ba0cceef8266bd627304ca`.
- [Archive review](run_result_v1/ARCHIVE_REVIEW.json) preserves the producer's
  complete row, guard flags, call counts and limitations.

This initial report transcribes saved results without scientific recomputation.
Completion acknowledgement and shared-status integration remain with root.
