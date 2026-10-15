# Flipped-basin relaxation scout — actual job 4831

The first declared source, A0ACD6B9F2 Ca-conditioned sample 0, passed the execution and physical-admission scout. Both metal searches converged, their final geometries passed the original local and full-source relative-contact checks, and both proposals were cross-scored under both metals. This supports root's decision about the remaining 18 declared searches. None were launched by this review.

## Actual response

All entries below are saved MACE vacuum works in kcal/mol relative to this source's original origin, copied from the completed collection.

| Geometry | Ca work | La work |
|---|---:|---:|
| Original origin | 0 | 0 |
| Original adaptive Ca | −8.516707 | −4.872677 |
| Original adaptive La | −7.206505 | −6.993732 |
| Fixed amide flip | +20.531480 | +43.031280 |
| Relaxed flip, Ca search | +11.913658 | +35.352249 |
| Relaxed flip, La search | +15.912722 | +28.707642 |

Common3, common4 and common6 therefore select the same original candidates under both mathematical and operational policies: adaptive Ca for Ca and adaptive La for La. The Ca-minus-La selected work is −1.5229743433870544 kcal/mol; both per-metal changes and the common4→6 paired change are exactly zero. The new searches lower their respective fixed-flip starts by 8.617822 and 14.323637 kcal/mol, but neither new geometry enters the selected pool minimum. The scout supplies no solvent-stage trigger or class-improvement result.

## Fidelity, admission and costs

- Ca: successful SLSQP termination, 13 iterations/14 function evaluations; maximum physical-heavy displacement 0.317553 Å; no boundary flag. La: 14/16, 0.406281 Å, no boundary flag. The exact original four modes/settings and cached flip-origin energies/forces were retained; gradients used the rebased Jacobian.
- Covalent/map/Jacobian checks passed. Both final full-source contact receipts pass with no newly severe contacts. Their general nearest-contact fields include covalent pairs; they do not establish the final separation or relief of the specific original short amide-H contact. Omitted-atom packing remains a limitation of the energy model.
- Actual execution: 569.195202 s on seven allocated CPUs, 3,984.366412 allocated CPU-s; 30 model forwards, zero failed worker requests, peak worker RSS 7,459,376 KiB. No new SCC, DFT or GPU work.
- The result retains all ten source rows. One source/two searches are complete; the other nine sources/18 searches remain unexecuted. All 80 prior endpoint receipts were checked and reused. The complete-ten-source solvent eligibility gate remains false.

These are bounded vacuum proposals. They do not establish an equilibrium rotamer, an unconstrained minimum, a composite stationary point, or predictive improvement.

## Immutable evidence

Workspace: `workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_scout_v1/`.

| Artifact | SHA256 |
|---|---|
| `CONFIG.json` | `053b7835217d8c51d24bab7e74e8781eff39acff8db24acb60b84d3e07845616` |
| `execution/FINAL_COLLECTION.json` | `e0d0e8e31eddf6954dfcf77ade9a871e871bdf67c8e688fb9f0f67fcc0ce758f` |
| `execution/EXECUTION.json` | `091ea56153952ddd0b1c978594fa2b1235085c94e0c06ac4b8ea2d621f994ce1` |
| `execution/mace/lane0/CPU_WORKER_SUMMARY.json` | `2812a4ee2864a99cb0388ddb7f4d408e77a86c294628eeac1bdb4b53a92e06e1` |

The collection pins both proposal receipts and prepared admission; the source result pins final full-source checks and every cross-cell. This report reads saved results only; it performs no new numerical analysis.
