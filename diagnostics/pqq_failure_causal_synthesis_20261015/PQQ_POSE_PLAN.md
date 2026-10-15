# Finite intact-PQQ pose response

Declared15October2026UTC after saved-force job4598 and before new geometries or
energies. This is a development test on consumed ExaA two-Ca sources0/3, not
validation or a claim that a source is wrong. Production stays unchanged.

## Motivation and fixed coordinate

Saved MACE forces show large shared quinone-shortening loads, which do not
specifically explain source discrimination. A separate source-dependent rigid
PQQ torque persists before and after donor accommodation. The current four donor
modes exclude PQQ entirely. Test this missing physical degree of freedom without
changing cofactor state, internal bonds or choosing a desired preference.

Use chemical e3 = e1 cross e2, with e1 along generic C5->O1 and e2 the normalized
component of C5->C6 perpendicular to e1. Center is the24-heavy-atom PQQ centroid.
Rotate the intact27-atom PQQ, including all three existing physical H, by exactly
-1 and +1 degrees. All other atoms remain fixed. No carbonyl shortening, ideal
bond target, new water, protonation, redox change, independent cap motion or
optimizer is included. One degree limits a several-Angstrom cofactor lever arm
to a small sub-Angstrom displacement; actual maxima must be reported.

Apply this SAME construction to each existing origin, adaptive_Ca and adaptive_La
candidate in both source contexts. Both endpoint metals score every admitted
geometry. This is24 new MACE vacuum endpoints plus12 exact existing baselines.
The complete shared pool is the original three plus their six rigid-pose variants.
Do not infer usefulness from a two-candidate subset or a favorable single row.

## Gates and calculations

Preparation is Slurm-only. Pin actual source maps, model and execution inputs.
Require identical Ca/La coordinates apart from metal identity, unchanged charge,
protons, inventory, cap mapping and atom order. Verify all PQQ internal distances
and existing H-parent distances are preserved within1e-10Angstrom; exact excluded
atoms/cap anchors must not move. Verify right-handed frame and deterministic
construction. Replay existing physical donor admission and contact rules. New
severe interfragment contacts use the inherited1.0Angstrom heavy/0.55Angstrom H
limits relative to the corresponding unrotated candidate, not a fitted cutoff.
Retain each failure; any missing required matrix cell leaves the full pool
comparison unavailable. Candidate counts are not populations.

Use the existing warm native MACE-OMOL float64 CPU executor and checkpoint,
with no new weights or dispersion terms. Exact twelve baselines are reused by
geometry/state/model identity. Record unrounded energies, force receipts, raw
per-metal work, geometry choices and common3->9 mathematical row-minimum shifts.
Report any inherited operational near-tie policy separately. No DFT or composite
forces are claimed; existing ExaA DFT validated donor work, not PQQ motion.

## Staged interpretation and cost

First run the finite MACE matrix. If neither source's full-pool paired contrast
changes by at least1kcal/mol in absolute value, close this finite-angle question
without adding solvent calculations. This is a declared practical resolution
for further work, not an uncertainty estimate or a label-based success gate.
Either direction qualifies; a worse result must remain visible.

If either source reaches that scale, evaluate the same24 new cells with the
qualified corrected-parameter standalone GFN2 and native frozen-vacuum-charge
ALPB recipe, with the exact existing numerical conventions. Reuse the12 old
solvent cells and compare complete composite pools. No native ORCA erroneous La
parameters, unavailable correction filled with zero, method replacement or hidden
baseline fallback. Record source changes and costs; no old-band class calls or
calibration adjustment. Any promising model must later face the designated
canonical references and full transfer population.

Use standard shared Slurm, the existing seven-CPU warmMACE layout, actual allocated
memory and numerical-thread bounds. No priority/memory/time override. This is
a fixed24-cell proposal test, not an open-ended pose search. Recent finite MACE
matrices took minutes; actual latency and allocated resources must be measured.
Root owns submissions and the stage decision. Curie owns the focused implementation,
Heisenberg reviews cofactor mapping/geometry, and Gauss remains on the independent
MMOL reference. No duplicate chemistry or production changes.
