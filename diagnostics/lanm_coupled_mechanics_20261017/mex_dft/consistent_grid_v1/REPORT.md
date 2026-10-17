# Mex consistent-grid qualification: execution checkpoint

Latest, 17 October 2026 UTC: central13281 is complete and admitted under the
original native rules. Read [the central report](CENTRAL_REPORT.md) for full
residual caveats, actual497700CPU-s cost and the approved two-singlepoint
continuation. Directional agreement still requires both displaced energies.
Read CURRENT.md and actual submission records for launch state. The startup
checkpoint below is historical; root owns execution.

## Why this test

The original Mex La external-charge check failed: analytic derivative
−3.673754 versus finite difference −4.087955 kcal/mol/Å. Saved diagnosis13277
attributes most of their −0.414201 difference to the final integration change.
The printed COSX contribution is −0.401606; subtracting it diagnostically still
leaves a residual outside the original tolerance. Raw point-charge gradient
mapping reproduces exactly. This suggests a numerical consistency test, not an
established physical model repair or ORCA implementation defect.

This separate recipe converges the density using the original final COSX grid.
It preserves the original geometry, electronic state, functional, basis, SCF
settings, finite-difference step and derivative tolerance. The exact two added
input lines are in [the execution specification](EXECUTION_SPEC_V1.md).

## Executed preparation

| Job | Result | Allocated cost |
|---|---|---|
|13279|Saved preflight failed on a header delimiter; zero molecular calls|10 CPU-s|
|13280|Saved preflight passed; zero molecular calls|10 CPU-s|
|13281|One fresh central La EnGrad; running|Pending actual accounting|

All three original outputs print four dots before the MaxIter value. The first
parser expected three. Versioned source `CONSISTENT_GRID_V2.py` changes only that
delimiter, retaining the exact125 requirement. Both source versions, independent
reviews and the failed result are preserved. The static review missed this
formatting mismatch; actual Slurm preflight caught it.

Preflight13280 replays the original failed collector exactly, validates archived
MM/ALPB reuse and all three original grid headers, and rejects each old grid as
the new candidate. Three prepared tasks preserve exact geometry/field/state
bytes. Only the declared grid lines are added. Root and independent review
checked the actual inputs. Result SHA256 begins `4c818c8d`, manifest `e0b4ce62`,
source `1600a0a5`; full pins are in [the archived result](preflight_v2/RESULT.json).

## Current execution and next gate

Central13281 uses config `28abed13`, one node and12 actual MPI ranks, with
scheduler memory/time defaults and allocation-derived MaxCore. Completion
watcher2759506 and alert-only health armer2759507 are recorded in
[the exact submission](central_submission/SUBMISSION.json). The native deck
has started. Actual startup review passes: grid2 and grid3 both use the original
final-grid construction, including1,026,059 points; XC/grid1/grid3 are unchanged.
Native charge0/multiplicity1/NEL832, MaxIter125, Total+1el-Energy/ConvForced0
and12MPI/MaxCore3100 also match. This validates startup settings, not a completed
SCF calculation or analytic gradient.

Both ±singlepoints require a completed, admitted central gradient and root's
review of its actual cost. At most three fresh native calls are allowed in total;
no native retries or new MM/ALPB calls. All three must use matching runtime
settings. The final comparison uses the fresh central derivative and both fresh
final energies under the original tolerance.

The original failed six-point continuation and dependent eight-cell response
panel remain stopped. A pass here would establish one La external-charge
direction at one step under this separate recipe; it would not establish
physical accuracy, electronic stability, global force accuracy or La/Dy affinity.
The independent whole-peptide reference13051 continues with its existing wakes.

## Preserved artifacts

- [Original source and static review](STATIC_REVIEW.json).
- [Narrow parser repair review](STATIC_REVIEW_V2.json).
- [Failed preflight and recovery decision](preflight_failure_v1/ROOT_REVIEW.json).
- [Actual passed preflight review](preflight_v2/ROOT_REVIEW.json).
- [Central declaration](central_submission/DECLARATION.json).
- [Actual native startup review and saved header](startup_review/STARTUP_REVIEW.json).

Source/preparation checkpoint is pushed in `15587c2`. Subsequent startup and
result records must be read before further submission; delayed notifications
never request repeated chemistry.
