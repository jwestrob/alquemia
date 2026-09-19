# Matched hydration mechanism test

Jacob approved the preceding proposal with “Yallah!”: compare an observed
water-present/water-absent state for both Ca and La, holding protein composition,
protonation and the native electronic method fixed. Reuse compatible archived
endpoints. The target is differential water–metal coupling, not a new affinity
classification or equilibrium occupancy calculation.

## Current source inventory

Both already consumed alpha-lactalbumin geometries have complete archived
repaired-v3 cores and endpoint pairs. 1F6S retains A:HOH211/212; 6IP9 retains
A:HOH310/322/326. Use every retained water in a separate single-water deletion,
without choosing waters from their scores. These are five matched squares from
one biological group, not five independent affinity observations.

The first integrity check found all ten water O–H bonds between 1.162862 and
1.189717 A. Source preparation used PDBFixer 1.12.0 / OpenMM 8.5.1 with no supplied
forcefield. The installed OpenMM fallback uses generic hydrogen bond springs
and repulsion, rather than a physical water forcefield. Parent heavy atoms are
unchanged; these hydrogen positions are preparation-generated, not crystallographic.

## Scientific choice pending

A question was sent to Jacob before new calculations: retain these archived
water geometries (10 new deletion endpoints, four archived full-state endpoints),
or first standardize all five water internal geometries using the exact existing
H2O-reference bond length/angle and preserve each O position, plane and bisector
(14 new endpoints: four full states plus ten deletions). Preserve all protein/cap
coordinates, charges, protonation and source oxygen positions in either path.
The recommended repair is not an orientation optimization or validation of a
water occupancy. Do not execute until this choice is recorded.

## Fixed calculations and interpretation

Native ORCA 6.1.1 r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 single points,
Ca charge -1 and La charge 0, both singlets. No functional/basis/ECP change,
SCF-tolerance search, geometry optimization, water replacement, multiple-water
deletions, new protonation state, reference refit or MACE inference. The common
water chemical potential and aquo offset cancel from each matched square.

DeltaS_add = [E_Ca(full)-E_Ca(deleted)] - [E_La(full)-E_La(deleted)].
Convert the final Hartree difference once with 627.509474 kcal/mol/Eh. Positive
means retaining that water shifts the descriptor toward La. Retain all endpoint
energies and available CPCM/SCF/D4/gCP components; do not assign a unique physical
cause from the CPCM term alone. Check direct-algebra closure to 1e-7 kcal/mol.
No absolute water-addition free energy, equilibrium population, calibrated
classification, new biological accuracy claim or entropy estimate is available.

Run through the existing affordable_workflow/run_orca_task_manifest machinery.
Recent matched source endpoints used 16 ranks and took 74–99 s each; the proposed
job uses 64 CPUs with four concurrent endpoints, expected scale several minutes.
This is an estimate, not a measured pilot cost. Preserve actual receipts and
failures. Jacob's no-compute-budget/no-time-limit direction remains in force;
record costs without an artificial stopping budget. No production/default change.
