# Donor accommodation survives native DFT in both PLM examples

**MACE identifies a real metal-dependent strain response.** In both compressed
PLM structures, moving the intact extra-Asp carboxylate away from its short
metal contact stabilizes La substantially more than Ca. Native r2SCAN-3c confirms
the direction and approximate magnitude. The known 4MAE control has a much smaller
response, also captured by the cheap models. All predeclared native calculations
are now complete; this is the final result of the fixed displacement experiment.

This gives the next classifier experiment a physical basis: account for how each
metal can accommodate nearby donor motion instead of allowing one compressed
source geometry to determine the contrast. It does **not** establish the PLM
proteins' true preferences, a relaxed minimum, or a corrected binding free energy.

## What worked, and what did not

Native OMOL reproduces all twelve endpoint-work signs and all six differential
signs. Its largest endpoint-work error is0.912kcal/mol on these selected paths.
The composite also preserves every sign, but its solvent contribution worsens
the differential agreement on both PLM examples. The largest composite endpoint
error is3.152kcal/mol. Thus the useful response comes primarily from MACE, not
from the extra solvent term. No general response-accuracy pass threshold was
declared; these observed errors must not be renamed a blanket qualification.

The solvent term remains part of the existing classifier and is reported intact.
Its differing behavior motivates testing whether MACE can propose donor
arrangements while solvent is evaluated at scoring endpoints, alongside the
already-running fully composite search. Neither route gets to choose geometries
by a desired class label.

## Complete native differential comparison

For each metal, work is E(q)-E(0) in the exact same expanded context. The
differential is work_Ca-work_La; a positive value means the displacement
stabilizes La more. The common reference offset cancels. Values below are
kcal/mol; [unrounded components](FINAL_COMPONENTS_v1.json) retain endpoint work,
native MACE, solvent transfer and composite separately.

| Case | Angle / rad | Native DFT | Native MACE | Composite |
|---|---:|---:|---:|---:|
| 4MAE | -0.2 | -1.476 | -1.686 | -1.613 |
| 4MAE | +0.2 | +1.087 | +1.014 | +0.948 |
| PLM8344 | -0.2 | -32.159 | -30.738 | -28.491 |
| PLM8344 | +0.2 | +29.868 | +28.535 | +26.957 |
| PLM07ab | -0.2 | -20.355 | -19.838 | -18.039 |
| PLM07ab | +0.2 | +20.335 | +18.917 | +17.674 |

Full PLM identifiers are PQQSEQ_83440678cbbd658047c9 and
PQQSEQ_07ab500e3df76b30d71c, original sample0 in each case. At the relief step,
the second PLM native Ca/La works are-12.314/-32.649kcal/mol. Its differential
error is-1.419 for MACE and-2.661 for the composite. The earlier first-PLM
component audit remains in [partial_v3](PARTIAL_DFT_REPORT_v3.md).

## Limits that determine the next experiment

The fixed tests preserve source scaffold outside the donor motion, PQQ, waters,
protonation, charge, multiplicity and physical cap mappings. They validate only
the chosen small extra-Asp displacements in three consumed contexts. There is
no native Glu/coupled-mode or relaxed-candidate validation in this experiment.

The cheap PLM profiles still descend to the+0.4rad boundary and their curvature
changes with interval size. Do not extrapolate a harmonic relaxation energy or
entropy from them. The separate bounded nonlinear pilot evaluates actual energies
and records interior, stationarity and curvature checks. A separate inexpensive
proposal experiment tests whether continuous solvent-gradient evaluations are
needed. Both retain all failures and original scores; neither changes production.

Both PLM proteins remain biologically unlabeled. All three saved source folds
of each were Ca-supported in the earlier sensitivity experiment, so choosing
another fold alone did not resolve the concern. A useful response must now earn
classifier utility under a common rule on established PQQ references. Existing
canonical fidelity and the ongoing matched DFT benchmark remain separate evidence.

## Execution and actual checks

All16new native endpoints completed normally with converged SCF; all12planned
displacement comparisons are available. Two compatible4MAE origin energies were
reused. Direct raw-output extraction matches all16 collected energies exactly.
No native retry, new geometry or changed method was needed. The original cheap
experiment retains63/64complete composite points: its two nativeGFN2 failures
at one1H4I La/Glu point remain missing, not repaired or omitted.

Native job1203771 used9214allocated wall seconds on64CPUs:
**589696allocated core-seconds**, zero GPU time. Batch MaxRSS was49525256KiB.
Together with the original cheap profiles, allocation cost was614944core-seconds
and13GPU-seconds. These are one-time development validation costs, not routine
scanner cost; local preparation/report/test CPU and old reused calculations are
outside that allocation total. [Receipts and scope](FINAL_EXECUTION_v1.json).

The immutable frozen analyzer produced
`workspaces/accommodation_torsion_20260920/final_DFT_result_v1.json` from the final
`prepared_v3/dft/collection_1203771.json`. It retains all64cheap points, seven paths,
exact replay controls and twelve native comparisons. No new fitted correction,
threshold, population, classifier accuracy claim or baseline change is present.

**Recommendation:** pursue physical accommodation using actual energy-selected
arrangements, preserve solvent-component scrutiny, and test the common rule on
known references. The mechanism is supported; discriminatory benefit remains
the next question.
