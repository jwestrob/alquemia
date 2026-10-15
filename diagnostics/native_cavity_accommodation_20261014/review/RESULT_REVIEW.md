# Independent result review — cavity-dependent accommodation

**Close this finite electrostatic-exclusion branch as a remedy for the ExaA
source discrepancy.** All existing common-pool selections remain unchanged.
Burial changes the differential accommodation only slightly; most of the
electrostatic augmentation is a shared origin shift. This conclusion concerns
these consumed sources, saved geometries and fixed-charge native model.

## Actual results

Values copied from root's Slurm 4587 saved-output summary; no numerical analysis
was executed by this reviewer.

| Source | Paired pool change, kcal/mol | Accommodation part, kcal/mol | Ca selection | La selection |
|---|---:|---:|---|---|
| sample0, control | +5.058230 | +0.026684 | adaptive_Ca, unchanged | adaptive_La, unchanged |
| sample3, target | +4.887095 | −0.049774 | adaptive_Ca, unchanged | adaptive_La, unchanged |

The saved target-minus-control change is −0.1711343629157227 kcal/mol;
−0.07645841008873511 kcal/mol comes from accommodation. Positive paired
Ca-minus-La shifts are La-directed. These diagnostic changes are not calibrated
class assignments or affinity changes. No threshold rescue or new geometry search
is supported by this result.

Native 4579 completed all 16 declared updates, with eight available endpoint rows.
Individual endpoint and per-metal work precision gates pass. Every paired and
source-level accommodation precision flag is false at the unchanged conservative
`1e-6 Eh` ceiling; the inherited origin source-contrast flag is also false. These
flags must accompany the reported sensitivities rather than be erased by matrix
completion. No precision rerun or relaxed ceiling is recommended.

## Comparator and compatibility review

The comparator uses the identical ordered common pool
`origin, adaptive_Ca, adaptive_La` for both metals. It checks the saved baseline
works and selections, adds only the declared electrostatic transfer, and separates
the origin shift from the change in origin-relative accommodation. The original
operational selections stay separate. Missing cells fail the summary; no partial
minimum or class band is introduced.

Independent metadata comparison verified **all 12** baseline
`evidence[metal][candidate].native_result` pins against the new/reused rows'
`archived_local` pins exactly. The baseline is midpoint 4137's immutable collection,
using only its `historical3` section; its midpoint/transplant candidates are not
included. The final comparator additionally verifies source IDs, metal/sample
identity and prepared-native pin joins. This was a metadata check, not a repeat
of the Slurm numerical summary.

The diagnostic still excludes full-source HB/SASA changes, permanent protein
charges and density response. Synthetic caps remain cavity centers; cap-volume
controls do not qualify a source-only physical cavity. Passing the local physical
and severe-overlap guards is not validation of a complete protein-solvent energy.
The negative result closes this exclusion mechanism at the tested scale, not all
possible environmental contributions.

## Measured execution and pins

The native-stage receipt records 205.32169128302485 wall seconds,
821.2867651320994 allocated CPU-seconds, 198.16 measured child CPU-seconds and
7,545,856 KiB peak child RSS. It contains 16 actual child execution receipts and
zero new SCC/MACE/DFT calls. Scheduler 206 seconds on 4 CPUs uses a different timing
boundary; neither figure is whole-campaign cost or a speedup measurement.

- Native `FINAL_COLLECTION.json`: `2fc07756da9c64857f012b957b4a01fd797e23ae6dbc364bcad5f2461b0239b2`.
- Native `EXECUTION.json`: `e4218d22f1e15d09353d458c926cbcf07ade3e8d73bccb087301ee994a06f6d2`.
- Summary `RESULT.json`: `87a86cac4135491361728ba172dceaf268acabd3f079343f5433ea75988aeaa0`.
- Exact baseline: `e7490cd21dd75ce81cde99660f48ecc541a2bffc6b39d5ba7aada2388116f6e7`.
- Reviewed `POOL_SUMMARY.py`: `1a1f7e67899357f77cbe370777740ffda531cfce832a8184904e77abf60c0260`.

Actual collections are under
`workspaces/native_cavity_accommodation_20261014/{native_v1,summary_v1}/execution/`.
The baseline is
`workspaces/exaa_metal_midpoint_20261014/scoring_v1/execution/FINAL_COLLECTION.json`.
