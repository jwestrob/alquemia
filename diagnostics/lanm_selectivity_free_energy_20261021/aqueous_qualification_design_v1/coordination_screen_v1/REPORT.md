# Nine-water transfer exposes a remaining water-model failure

**The frozen water candidate fails three of four declared checks.** The two
calculations completed successfully; this is a scientific transfer failure.

| Check | Model/reference comparison | Existing limit | Result |
|---|---|---|---|
| Paired eight-to-nine-water energy change | Model −5.485; reference −0.534; error −4.951 kcal/mol | Absolute error ≤1 | Fail |
| La rigid response | RMS 3.811; maximum component 7.961 kcal/mol/Å | RMS ≤2; maximum ≤5 | Fail |
| Dy rigid response | RMS 4.686; maximum component 9.949 kcal/mol/Å | RMS ≤2; maximum ≤5 | Fail |
| Dy-minus-La rigid response | RMS 1.054; maximum component 2.268 kcal/mol/Å | RMS ≤2; maximum ≤5 | Pass |

The paired energy is `[Dy−La] at W9 − [Dy−La] at W8`. Free-water energy
cancels across the metals. The original W9 geometry also rearranges the original
eight waters; this is a coordination-and-geometry transfer test, not an isolated
water insertion or a prediction of equilibrium hydration populations.

The reference is the existing provisional frozen-f ECP-PBE0-D4 surface. Its
remaining physical-state and method uncertainties are retained. This discrepancy
alone does not identify a unique missing interaction or establish DFT as truth.

## What changes

The previous fitting controls and four independent orientation results remain
valid for their tested configurations. They did not establish coordination
transfer. The BIO namespace check also remains passed; it preserves the intended
model, including the explicit water internal-term difference.

Bulk sampling with this candidate is not released. A saved comparison with the
original model is being prepared to determine whether this discrepancy predates
the water repair and which model terms contribute. No parameter values, gates or
live donor-fit settings have changed; no refit follows automatically.

## Execution and preservation

Job 17670 completed 0:0: two new fixed model cells, four native programs,
zero new DFT or finite-difference calculations, three allocated CPU-seconds.
Both W8 model anchors and all four quantum references were reused.
All inputs, raw outputs, comparison and root review are in the corresponding
workspace `coordination_screen_v1/run_v1`. Independent actual review is pending.

## Documentation defect

Root accidentally overwrote the pre-submission review file with the terminal
review. The original review pin in the preserved release therefore no longer
matches. This defect is recorded in `run_v1/PROVENANCE_NOTICE.json`; the lost
review has not been reconstructed. Source, bound inputs, release/submission
records, raw native outputs and scientific result remain intact. Independent
scientific review must retain this provenance limitation. No chemistry is rerun
or historical release rewritten to conceal it.
