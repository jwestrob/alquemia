# Hydration-square result — 2026-09-18

**All 14 ORCA endpoints completed. Water preparation and identity both materially
affect the electronic Ca−La contrast. This is a hydration-mechanism result, not
a demonstrated accuracy gain or a released occupancy correction.** Baseline
inputs, archived results and production defaults are unchanged.

## What changed and why

The approved pilot uses two already-consumed bovine alpha-lactalbumin structures,
1F6S and 6IP9, one biological group. All five observed waters were included,
without selecting favorable outcomes. The generic preparation had stretched
water O−H lengths (1.163–1.190 A). In a new path we normalized each water to the
pinned H2O reference while preserving its O, plane and bisector and every other
atom coordinate. This normalization is not an orientation optimization.

Both full parents and every individual water deletion were evaluated for Ca
and La: seven pairs, 14 native ORCA 6.1.1 r2SCAN-3c/CPCM(Water)/DefGrid3 single
points, same archived SCF policy, Ca −1/La 0 singlets. No MACE call, optimization,
added water, altered protonation, new aquo reference or threshold fitting.
Protocol: `native_r2scan3c_cpcm_single_water_square_v1`, preparation policy
`reference_internal_geometry`. See [decision](DECISION.md), [source finding](PREPARATION_FINDING.md)
and [execution configuration](EXECUTION_CONFIG.json).

## Separate effect of repairing water hydrogens

Full-parent R=E_Ca−E_La changes by **+4.925351 kcal/mol (1F6S)** and
**+2.479708 kcal/mol (6IP9)** relative to the archived stretched-water parents.
Positive is more La-like. Thus the preparation defect affects the contrast,
even with identical protein and oxygen coordinates. It does not by itself
establish a correct water orientation, occupancy or a repaired classification.

## Effect of each water after the common repair

`DeltaS_add = [E_Ca(full)−E_Ca(deleted)] − [E_La(full)−E_La(deleted)]`.

The same water chemical potential and aquo offset cancel. Convert the final
Hartree contrast once using 627.509474 kcal/mol/Eh. Positive means retaining
that water shifts the contrast toward La; negative means toward Ca.

| Structure | Water | Metal−O distance (A) | DeltaS_add (kcal/mol) |
|---|---|---:|---:|
| 1F6S | A:211 | 2.577 | -19.070 |
| 1F6S | A:212 | 2.581 | +3.683 |
| 6IP9 | A:310 | 2.556 | -19.951 |
| 6IP9 | A:322 | 2.702 | +5.008 |
| 6IP9 | A:326 | 2.483 | -0.452 |

Retaining A211 in 1F6S or A310 in 6IP9 shifts the contrast about 19–20 kcal/mol
toward Ca. A212/A322 instead shift it 3.7–5.0 kcal/mol toward La; A326 has a
smaller effect. All five outcomes are retained. These are conditional effects
with all other waters present, **not additive contributions** or five independent
biological validation cases. Atom identities and the metal−O distances above
come from the pinned parent mapping/coordinates.

This supports explicitly representing water identity and local geometry rather
than using water count alone. Removing a water also changes continuum cavity and
electronic density, so the effect is not uniquely an isolated metal−water bond.

## Energy components

Each entry is the same Ca−La matched water-addition contrast in kcal/mol.
SCF **already includes CPCM**. Total is SCF+D4+gCP; adding CPCM again would count
it twice. The printed CPCM term cannot uniquely isolate a cavity or polarization
mechanism from coupled electronic changes.

| Structure/water | SCF | CPCM (included in SCF) | D4 | gCP |
|---|---:|---:|---:|---:|
| 1F6S A:211 | -18.510 | -0.530 | +0.032 | -0.592 |
| 1F6S A:212 | +3.791 | -1.769 | +0.020 | -0.128 |
| 6IP9 A:310 | -19.320 | +0.629 | +0.029 | -0.661 |
| 6IP9 A:322 | +5.084 | -1.940 | +0.017 | -0.093 |
| 6IP9 A:326 | -0.180 | -2.276 | +0.019 | -0.291 |

## Actual energies and validation

Unrounded parser values, component terms, all endpoint hashes and receipts are
in the [result index](RESULT.json) and its pinned collection. Values below retain
the 12 decimals printed by ORCA; units are Hartree.

| State | E_Ca | E_La |
|---|---:|---:|
| 1F6S__full | -1855.780459462680 | -1209.692433436387 |
| 1F6S__minus_A_211 | -1779.390824319824 | -1133.333187487176 |
| 1F6S__minus_A_212 | -1779.338320259036 | -1133.244425635715 |
| 6IP9__full | -1932.160660232464 | -1286.063982380329 |
| 6IP9__minus_A_310 | -1855.777886841956 | -1209.713003689456 |
| 6IP9__minus_A_322 | -1855.737201409227 | -1209.632542288746 |
| 6IP9__minus_A_326 | -1855.739396974758 | -1209.643439995769 |

Four prior real-fixture software tests passed (0.341 s); all seven prepared pairs
passed exact paired-coordinate, retained-atom, neutral-water-deletion, charge,
parity and reference-geometry checks. The existing runner dry-run passed. All 14
actual electronic calculations converged and terminated normally, no retries.
All five square closures pass the predeclared 1e-7 kcal/mol algebra tolerance.
Original archived parents and energies were reverified after execution.

This algebra tolerance is not an electronic-structure accuracy estimate. No
SCF/grid refinement or experimental occupancy validation was performed. The
existing production PQQ calibration was not rerun or modified by this pilot.

## Measured cost

Job 1201801, node-64-768g-6, four 16-rank endpoints concurrently, 64 allocated CPUs,
256 GB requested, no GPU. Scheduler: 333 s wall, 21,312 allocated core-seconds,
reported TotalCPU 04:45:45 and batch MaxRSS 7,525,768 KiB. Runner portion:
331.638098605 s wall / 21224.838310719 allocated
core-seconds. Endpoint wall range 68.617–109.820 s, median
88.143 s. Per-endpoint assigned-rank sum 19740.381328 s.
Scheduler MaxRSS is the reported batch accounting field, not a measured global
peak across independently monitored processes. No compute-budget stop or
user-requested time limit; scheduler defaults applied. Earlier audit/test costs
are recorded separately; no earlier quantum attempt belongs to this pilot.

## What this enables next

The useful next modeling target is metal-dependent hydration-state selection
with a consistent bulk-water reference and a controlled treatment of water
orientation/relaxation. These squares do not determine whether either metal
actually expels or retains a given water: individual addition free energies do
not cancel the water reference. Do not choose a favorable deletion because it
matches the assay, sum these conditional deletions, or insert them as a blanket
classifier correction. First define the complete exchange cycle and state
selection rule, then evaluate both metals on the same declared candidate states.
No extra state calculations were launched by this report.

**Recommendation:** pursue explicit hydration-state modeling; retain the working
baseline while it is developed. This inexpensive pilot identifies a consequential
preparation defect and a substantial water-specific coupling worth modeling.
It does not yet show improved predictive accuracy. Replay/collection commands
are in [COMMANDS.md](COMMANDS.md). All jobs from this pilot are finished.
