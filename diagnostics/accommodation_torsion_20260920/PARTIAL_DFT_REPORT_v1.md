# Partial native DFT adjudication: 4MAE control complete

**The four completed control displacement energies agree with the composite on
all work signs; differential-work errors are 0.137–0.139 kcal/mol.** This is useful
local evidence for the small reference motion. The large predicted relief in the
two compressed PLM structures remains untested until their native energies finish.
No new scientific calls, changed angles, thresholds or execution interventions
were introduced for this snapshot.

| Extra-Asp angle | DFT Ca work | Composite Ca | DFT La work | Composite La | DFT Ca−La work | Composite Ca−La |
|---|---:|---:|---:|---:|---:|---:|
| −0.2 rad | −0.6225 | −0.3086 | +0.8536 | +1.3048 | −1.4762 | −1.6134 |
| +0.2 rad | +1.6223 | +1.8935 | +0.5349 | +0.9452 | +1.0874 | +0.9483 |

All values are kcal/mol relative to each metal's exact reused 4MAE context
origin. Maximum individual endpoint-work error is 0.4512 kcal/mol. The error in
Ca-minus-La work is −0.1373 at −0.2 rad and −0.1392 at +0.2. No fitted correction
was applied. These values support this particular control response; they do not
establish broad force accuracy, a relaxed minimum, entropy, or classifier gain.

## Execution and missing values

Job 1203771 remains active. At this snapshot:

- 4/16 new native endpoints complete, all with normal termination and SCF
  convergence verified through their actual execution receipts.
- Four PLM endpoints have active outputs and no observed SCF-failure or error-
  termination marker; eight have not started. All 12 retain unavailable native
  energies in the immutable partial collection and analysis.
- The first completed endpoints took 26.7–45.0 minutes each on their 16-rank task
  allocations. The existing four-task concurrency is operating normally; these
  are expanded-context development calculations, not proposed routine scanner
  costs. Final allocation cost is still pending.
- The earlier two GFN2 failures at the unrelated 1H4I La Glu −0.2 point remain
  visible. No failed energy was filled or retried. No DFT execution issue needs
  attention based on this snapshot.

## Immutable artifacts

- `workspaces/accommodation_torsion_20260920/prepared_v3/dft/partial_collection_1203771_v1.json`:
  four complete native values, full 16-task denominator, pending values null.
- `workspaces/accommodation_torsion_20260920/partial_DFT_result_v1.json`:
  unchanged frozen analyzer applied to the partial collection; four available
  native work comparisons out of the 12 preselected displacement comparisons.
- [PARTIAL_EXECUTION_v1.json](PARTIAL_EXECUTION_v1.json): actual endpoint receipt
  and running/not-started status snapshot.

The original plan, primary report, production scorer and all prior results remain
unchanged. The next action is final collection and PLM adjudication after the
existing 16-task job finishes; no additional calculation is queued here.
