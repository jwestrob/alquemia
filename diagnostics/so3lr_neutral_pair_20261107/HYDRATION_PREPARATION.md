# Bounded hydration preparation repair

Declared 7 November 2026, after saved-data job53123 and before this repair.
All six original endpoints have zero water oxygens within3.5 A of either metal.
The ligand contacts differ, but the original three boxes do not provide the
intended diversity of initial hydration. Preserve them and every earlier
potential result. The queued periodic benchmark52844 keeps its original inputs.

Prepare one new set for the subsequent dynamics scout. Keep all original
ligand conformations, stereochemistry, metal positions, counterions, atom order,
charge, spin, cell and total water count. Insert two waters per metal for
starts0/2 and three for start1, then repack the remaining bulk water. These are
initial solvent arrangements, not constraints on later coordination or denticity.

Use the exact water geometries and metal-relative positions from archived
aquo8_central_La for starts0/2 and aquo8_central_Dy for start1. Both A/B endpoints
share the same coordinates, including their hydration. For each site, enumerate
the24 proper cube rotations combined with six global z rotations(0,15,...,75
degrees), and every subset of the required water count from the eight source
waters. Choose the subset/orientation with the largest minimum atom clearance
from the existing solute and its other selected waters. Require at least1.6 A.
This geometric selection does not use model energies, experimental preferences
or selectivity residuals. Preserve the complete chosen transformations and source
hashes. No new ligand optimization is performed.

Fix those starting waters only during Packmol preparation; use the pinned
Packmol20.15.0, periodic42 A cell,2 A packing tolerance and seeds55100–55102.
Pack2396/2394/2396 additional waters, respectively, for2400 total. Restore only
file-rounding changes to the fixed coordinates. Require no atom overlaps below
1.6 A between independent fragments, finite unchanged ligand geometries and
exact endpoint correspondence. Re-run the already declared initial observations
on the repaired boxes. Do not claim the resulting starts are equilibrated.

Allocate three CPUs with default memory/time under standard_expand. An internal
180-second timeout plus five-second grace bounds this single preparation
attempt and its initial observations to555 allocated CPU seconds. No automatic
retry or further conformer/hydration search follows a failure. This remains
within the existing four-CPU-hour CPU-only preparation/static allocation; the
separate one-H200-day envelope includes its necessary host CPUs and all GPU
benchmarks/failures. Record both resource types without hiding host CPU cost.

Before MD on the repaired input, evaluate its start0 A/B, a shared rigid
translation/wrap, and one shared remote-water deletion: exactly six new L
evaluations with the original settings and 0.5-kcal/mol robustness criterion.
`periodic.py --hydrated-scout` selects this explicit branch. Its standalone
launcher has a360-second internal bound and is not submitted before review of
52844. That existing job's default execution path and inputs are unchanged.
Count these calls and necessary host CPUs within the one-H200-day allocation.
