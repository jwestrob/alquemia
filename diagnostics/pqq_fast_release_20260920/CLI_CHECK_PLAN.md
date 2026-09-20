# Literal standard CLI integration check

Recorded before execution, within Jacob's source-to-score promotion request and
explicitly approved by parent on 2026-09-20. After the full 28-case release gate
passes, run one additional fresh 1H4I through the actual
`affordable_workflow.py standard` entrypoint. This checks dispatch/orchestration
added after the full-panel harness was frozen; the harness alone cannot establish
that the public entrypoint works.

Exactly one fresh raw-source normalization/protonation/context preparation,
2 native MACE scalar calls and 4 native GFN2 calls. Same source, selectors, physical
state, checkpoint, numerical settings, composite expression and frozen bands.
No archived preparation or energy is required in the example request. Compare
coordinates and score with the already consumed 1H4I result after execution.
No extra DFT calculation, parameter fit, new labels or scientific variants.

Use the same working GPU allocation: one H200, 32 MPI slots/CPUs and 200000 MiB.
Retain actual source-to-score and allocation receipts. Failure blocks final
standard-mode publication until a technical recovery of the same inputs passes;
no scientific fallback or threshold padding. Additional runs count separately
from the full 28-case check: release totals are 29 fresh preparations, 58 MACE
and 116 GFN2 calls if this one check succeeds on the first attempt.
