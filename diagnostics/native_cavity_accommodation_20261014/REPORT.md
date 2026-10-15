# Cavity shielding does not improve this ExaA accommodation comparison

The completed finite test leaves all four selected source/metal geometries
unchanged. Cavity-dependent accommodation changes the source contrast by only
−0.076458 kcal/mol. Including the origin shift gives−0.171134 kcal/mol, primarily
from a common approximately+5 kcal/mol metal contrast. No classification gain
is established. **Close this electrostatic-exclusion branch; do not promote it
or recalibrate its common shift.**

Both sources and every declared cell are retained. All16 new updates and all
per-endpoint gates pass, with origins and local energies reused exactly. The
aggregate charge-rounding precision flags remain false; no tolerance was changed.
These are raw conditional sensitivities, not a newly qualified total solvent
model or a binding-affinity result.


All values below describe the declared electrostatic-only diagnostic. They
are not new class assignments or a complete protein-solvent correction.

| Source | Paired pool change | Accommodation part | Ca selection | La selection |
|---|---:|---:|---|---|
| sample0 | +5.058230 | +0.026684 | adaptive_Ca → adaptive_Ca | adaptive_La → adaptive_La |
| sample3 | +4.887095 | -0.049774 | adaptive_Ca → adaptive_Ca | adaptive_La → adaptive_La |

Energy changes are kcal/mol; positive Ca-minus-La changes are La-directed.
Source3-minus-source0 paired change: -0.171134 kcal/mol.
Accommodation contribution to that change: -0.076458 kcal/mol.

Conservative aggregate precision flags are retained in RESULT.json and source
collections. Passing molecular execution/per-endpoint gates does not make those
aggregate flags pass. Native HB/nonpolar burial and permanent environment charges
remain excluded. Synthetic-cap controls do not establish a source-only cavity.

No new electronic endpoints, geometries, thresholds or production changes.

## Evidence, cost and boundaries

Actual preparation4577:8/8,14s/1CPU,10.262 measured TotalCPU seconds,215884KiB
peak RSS. Native4579:16/16,206 scheduler seconds/4sharedCPUs (824allocatedCPU-s),
203.781 measured TotalCPU seconds,7662672KiB scheduler peak RSS. Its runner
boundary records205.321691s,198.16 childCPU-s and7545856KiB child peak; these
measure different intervals. Saved comparison4587:1s/1CPU,0.190CPU-s,29900KiB.
There were no new SCC, MACE, DFT, optimizations, folds or GPU evaluations.
The qualified native engine is serial; four slots provided default memory.

The original12-update scout4562 is separately preserved in
../native_cavity_shielding_20261014/REPORT.md. Existing local energies and charges
were used without conversion or state changes before native replay. Preparations
preserve all source atoms, mapping, paired coordinates, cap anchors and saved
physical donor motion; newly introduced severe full-source contacts are checked
using the inherited1.0Å heavy/0.55Å H policy, not a fitted new cutoff.

The independent review verifies all12 old baseline cell pins against actual
common-pool evidence. See review/RESULT_REVIEW.md and RESULT.json. Native full
components, cap controls, uncertainty bounds, partial receipts and exact library
pins are linked from EXECUTION.json and the workspace collections.

No statement here excludes permanent protein electrostatics, responsive electronic
populations, omitted hydrogen-bond/nonpolar burial, or different physical states.
Those are different models. This finite excluded-volume diagnostic does not earn
routine cost on the two consumed sources. No automatic wider source or parameter
sweep follows; PQQ production/defaults and previous failures remain unchanged.
