# 7CCN whole-peptide reference: scout admitted, three cells running

19 October 2026 UTC. Scout14795 completed the declared before-La native PBE0
energy and 256-atom gradient under the original admission rules. Root reviewed
the actual output and scheduler cost; independent review found no continuation
blocker. The remaining three fixed cells are running:

| Cell | Slurm job | Allocation |
| --- | --- | --- |
| after-La | 14880 | 12 MPI ranks, 48,000 MiB |
| before-Dy | 14881 | 12 MPI ranks, 48,000 MiB |
| after-Dy | 14882 | 12 MPI ranks, 48,000 MiB |

These complete the original four-cell panel. They reuse the frozen geometry
pair, method, parameters and solvent components. Each has one native invocation,
allocation-derived MaxCore3100MB/rank, and completion plus alert-only health
wakes. No memory/time override or waiting collector allocation was requested.
All three passed their startup preflight/scout replay and reached native launch.

## What the scout establishes

The native calculation reports normal termination and energy-check convergence
after15 cycles, with the complete finite gradient and original state/receipt
checks. Its result is `5fe74c9b`, native receipt `01571ff9`, output `457c1c6d`,
and gradient `685cc58c`. One native call was made; no new MACE, solvent, GFN,
MM or optimization calls were made by this panel so far.

The saved initial La force comparison gives MACE-minus-PBE0 atom-vector RMS
3.513854 and maximum8.128667 kcal/mol/Angstrom. Both initial composite force
projections point along the frozen deformation chord: PBE0+ALPB40.285425 and
MACE+ALPB38.510768 kcal/mol/Angstrom. These describe one nonstationary endpoint;
they do not establish the deformation energy or the La/Dy difference. All four
cells are required for that comparison, with the existing1kcal/mol reporting
resolution unchanged.

## Limitations retained

Only the energy residual meets its printed tolerance. Density, orbital-gradient,
orbital-rotation and DIIS residual flags remain false, with a small-gap warning.
This satisfies the declared native convergence semantics; it does not establish
physical-state stability or independent force accuracy. Those qualification
flags remain false. Vacuum-dianion, frozen-f Dy and post-added-solvent limitations
are unchanged. Missing Dy MACE forces remain unavailable.

7CCN is Lu-conditioned and uses the same peptide as7CCO. The existing80-step
terminal is nonconverged. This is a distinct fixed-deformation test, not an
independent-protein validation, affinity prediction or a new preference result.

## Cost and continuation decision

Scout14795 used23,440 scheduler seconds (6h30m40s) on12 CPUs:
281,280 allocated CPU-seconds. This lies within the prior whole256 reference
range23,058–24,919s. Native/application timing and measured scheduler CPU time
are preserved separately; allocated cost is not utilization.

Root's immutable decision `74c11c68` approves only after-La, before-Dy and
after-Dy, using actual scout and original preflight pins. Each task replays those
gates before execution. The first scout is never repeated. Collect the saved
results only after their completion; failed or missing cells prohibit a survivor
work comparison. The broader mechanics goal remains unmet.

Evidence: [independent review](SCOUT_RESULT_REVIEW.md),
[decision](scout_review_v1/CONTINUATION_DECISION.json),
[cost](scout_review_v1/SCHEDULER_COST.json),
[lossless archive pins](scout_review_v1/EVIDENCE_PINS.json),
[submissions](continuation_v1/SUBMISSIONS.json).
