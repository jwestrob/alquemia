# 7CCN-La: model force convergence achieved

Run15525 completed with **`termination=full_force_converged`**. The unchanged
whole-MACE + fixed-FF ALPB model met both full-force criteria after
**553 additional accepted steps and 567 additional distinct evaluations**,
reaching **725 accepted steps / 739 evaluations** cumulatively.

The final maximum atomic force is **0.07505475 kcal/mol/Å** and the RMS
atomic-vector force is **0.02961095 kcal/mol/Å**. Both satisfy the existing
thresholds of 0.1 and 0.03, respectively. All saved trace and trial geometry
guards pass, and the original source was not rebased.

**This establishes convergence under the model's force criteria. It does not
establish physical accuracy or reference validation of this segment.**

## Saved model results

| Full256 composite force summary | Step172 | Step725 | Unchanged threshold |
|---|---:|---:|---:|
| Maximum atomic force, kcal/mol/Å | 4.824511 | 0.075055 | 0.1 |
| RMS atomic-vector force, kcal/mol/Å | 1.045572 | 0.029611 | 0.03 |

The final full-force L2 norm is **0.47377515 kcal/mol/Å**. The termination was
force convergence, rather than exhaustion of the declared step/evaluation
budget. These endpoint summaries do not imply that forces decreased monotonically.

| Incremental work, step725 minus step172 | kcal/mol |
|---|---:|
| Whole-MACE vacuum | +11.053362 |
| Fixed-FF ALPB | −24.547593 |
| Composite | −13.494231 |

Within this model, solvent lowering outweighs the increase in vacuum energy.
This component accounting does not establish the accuracy of either term.
The producer reports cumulative work from the original start of
**−170.954464 kcal/mol** for the composite; this is a model result, not an
independent reference measurement.

## Scope and remaining limits

This was one La search on the consumed, Lu-conditioned 7CCN peptide, with all
256 physical atoms movable. Original bond, chirality, peptide-basin,
severe-overlap and represented-state guards were retained. The run made
**567 MACE calls and 567 ALPB calls**. It made no DFT, GFN, MM, opposite-metal
or other-search calls.

The earlier independent reference concerned step80→92. It does not validate
step172→725 or the final geometry. Physical-state qualification, force-accuracy
qualification and the scientific accuracy gate remain false. Conditional
terminal chemistry, the vacuum electronic model with post-added fixed-FF
solvent and the existing reference-state/SCF limitations remain relevant.

Small model forces do not prove a stable minimum, physical equilibrium,
accurate forces, binding affinity or metal preference. Dy force availability
is unchanged. The producer retains `no_automatic_followup=true`; root owns
the next scientific decision.

## Actual allocation and archive

| Job | Elapsed seconds | CPUs | GPUs | Allocated CPU-seconds | Allocated GPU-seconds |
|---|---:|---:|---:|---:|---:|
| 15525 | 406 | 12 | 1 H200 | **4,872** | **406** |

These are scheduler allocation costs, not utilization. Prior searches,
preflights, audits and references retain their separate costs. Application
elapsed time was **403.841720 seconds**; the producer records 6.410947 seconds
model loading, 69.992559 seconds in MACE calls and 28.555505 seconds in ALPB
calls. [Exact accounting](run_result_v1/COST.json) retains the scheduler record.

The [archive manifest](run_result_v1/EVIDENCE_PINS.json) preserves the complete
execution tree: **567 evaluation directories**, all raw component receipts and
gradients, **555 checkpoint files**, the complete **554-row trace**, saved parent
admissions and exact source/configuration/launch evidence. A single lossless
tar.xz stream avoids duplicating cumulative checkpoint content across thousands
of separately compressed files. Every member has an original byte count and
SHA256, verified by decompression; workspace originals remain unchanged.

- Actual RESULT SHA256:
  `95d870ebbfdedede5166c7499e4b8181dd8c24289503bdb2e35eaba02c4a5068`.
- Actual SEARCH SHA256:
  `a9495fd5f42c6fc115a26e00cd33dd9d2369ce4b7ffcdcb4fe5aa6ee0115c25d`.
- Final POINT SHA256:
  `6605f43cc3289509cee65c9be97a195abaad43d2d247ab488f5de889aa6f5419`.
- [Archive review](run_result_v1/ARCHIVE_REVIEW.json) retains complete producer
  flags, guard status, parent lineage and scientific caveats.

This report transcribes saved results without scientific recomputation.
No submission, source change or completion acknowledgement was made here.
