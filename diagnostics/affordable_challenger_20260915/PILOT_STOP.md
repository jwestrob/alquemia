# Jobs finished; environmental validation incomplete

Updates PROGRESS_RECOVERY.md. No jobs remain running for this pilot.

- Six successful quantum endpoints and six passing ESP-quality checks.
- APBS completed one six-solve state, qm33 La. The remaining 17 states were
  explicitly unrun after budget exhaustion. No paired environmental score,
  identity, refinement, transformation or partition result is available.
- Its coarse-grid homogeneous Coulomb cross check differs from the analytic
  value by **−1.3365592017929906 kcal/mol**, exceeding the frozen 0.5 tolerance.
  This is a failed coarse-grid check; convergence has not been tested. It is
  not a predictive classification result or a reason to change the tolerance.
- Full solver-stage accounting, including the parser-failed first attempt:
  **57,792 allocated core-seconds versus the 52,976 budget**. Aggregate across
  all three jobs: 153,216 allocated core-seconds. The in-flight solve overran
  its admission estimate; further solves were not launched.
- **Execution mistake:** one serial APBS process used 143.82 user + 11.13
  system CPU seconds over 155.9623 wall seconds, but occupied a full 344-CPU
  allocation. This poor utilization inflated allocated cost. The experiment
  does not establish intrinsic model unaffordability. A continuation needs
  measured-cost scheduling that batches independent solves within an allowed
  allocation, retaining all prior costs and the frozen scientific checks.

No additional calculations were launched. Fixing execution utilization is the
next engineering step before proposing a costed continuation. The baseline
remains default because validation is incomplete, not because the challenger
has been scientifically rejected. No new DFT endpoints are needed to continue
these checks; preserved wavefunctions and validated charge records exist.

Exact receipts and status: [PILOT_STOP_RESULT.json](PILOT_STOP_RESULT.json).
