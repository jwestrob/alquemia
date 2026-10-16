# Fixed-source, matched one-water LBT3 preparation

Declared 16 October 2026 before placement calculations or new energies.
Protocol: `lbt3_fixed_source_one_water_geometry_v1`. Root owns execution and
the later energy plan. This declaration authorizes geometry preparation only.

## Question and chemical inventory

Prepare a finite test of explicit water competition near Asn5, absent from the
current dry explicit inventory. Use BOTH original admitted 256-atom sources
from preparation12469, 7CCO and7CCN. Freeze every original coordinate, including
target metal, deposited H and completed termini. Add exactly one neutral H2O
using the actual archived geometry
`lanthanide_series_capability/calculations/constants/xyz/water.xyz`.
Every pose has259 atoms, charge−2, target index255, La multiplicity1 or Dy6;
the new O/H/H indices are256/257/258. Preserve original source graphs and maps.
No cap, proton change, counterion, peptide motion or missing-H minimization.

The water position is MODELED, not deposited or experimentally determined.
The original NMR structures used identical2.1–3.0Å Ln–O restraint bounds and
ion-specific refinement (archived primary paper Par24–25). Par12–15 discusses
water/Asn5 competition; Supplementary Fig.S10 locates a water RDF peak at2.6Å.
The currently declared terminal/protonation hypothesis remains conditional.

## Deterministic placement, frozen before geometry results

Use the exact mapped atoms `metal`, `A/5//OD1`, `A/5//CG` in each source.
Let e0 be the unit metal→OD1 vector. Let e1 be the normalized projection of
OD1→CG perpendicular to e0, and e2=e0×e1. A degenerate frame is unavailable;
there is no fallback atom choice or source alignment that moves the peptide.

Enumerate24 rays in order: polar angle60°, then75° from e0, with azimuth
0°,30°,...,330° around e0 in each ring. Direction is
cos(theta)e0+sin(theta)[cos(phi)e1+sin(phi)e2]. This probes the Asn5 side of the
site away from the occupied OD1 direction; it does not prescribe displacement
or decoordination of Asn5.

For EVERY ray construct the same four poses in this order:
inner/plane0, inner/plane90, outer/plane0, outer/plane90.
Water oxygen lies2.6Å (inner) or4.5Å (outer) from the fixed metal. H-bisector
points radially away from the metal. In plane0 the H-difference direction is
the normalized projection of e0 perpendicular to the ray; plane90 rotates
that direction90° around the ray. Apply a proper rigid rotation and translation
to the pinned real H2O; preserve its complete O–H–H geometry.

Select ONE ray per source by maximizing the minimum normalized new-water
contact clearance across ALL four poses and all original256 atoms. The
normalized clearance is distance divided by its fixed cutoff below. The metal
O distance is prescribed and contributes no selection score. Choose the first
candidate in the declared order on an exact floating-point tie. This is a
bounded geometric placement rule, not energy/force optimization. Save all24
candidate clearances and individual pass/fail results for review.

## Admission and no replacement

New-water/original-atom minimum distances, Å:

| Pair | Cutoff |
|---|---:|
| water O / peptide heavy atom |2.0|
| water O / peptide H |1.2|
| water H / peptide heavy atom |1.2|
| water H / peptide H |0.8|
| water H / target metal |2.0|

The2.0Å O/heavy cutoff reuses the existing `hydration_site_proposals` flag;
the H limits are predeclared conservative new-placement screens, not fitted
energetic radii. Record nearest contacts and every violation. Apply these
screens only to new-water contacts: inherited peptide contacts are unchanged.
Require exact original-coordinate equality, paired-metal equality except the
target symbol, complete neutral-water inventory, correct state/electron parity,
proper rotation and all three water pair distances preserved within1e−10Å.
Allow1e−10Å numerical roundoff at declared contact/radial cutoffs.

Keep the four slots belonging to the selected ray even if some/all fail.
Write and label rejected geometries; do NOT replace a failed orientation,
change the ray after admission, alter a cutoff or insert another water.
No admissible inner/outer pair is an explicit limitation of this fixed-source
test. Preparation does not prove physical hydration occupancy or equilibrium.

## Artifacts and later interpretation

Emit exact config/implementation/plan/template/source/helper pins, source-frame
atom mapping and vectors,24-ray audit, four physical pose/map records per source,
and16 flat paired-metal task slots with explicit `admitted` flags. Whole-source
and water bond maps remain separate; no metal–water covalent edge is invented.
Reuse `hydration_basin_coordinates.WaterCoordinates` for rigid transforms.
Use current MACE/solvent physics only through a separately reviewed259-atom
adapter after preparation admission; no old256-atom energy receipt is reusable
for the changed inventory. This preparation makes zero model/backend calls.

The later finite panel is at most16 composites (16MACE,32GFN2 calls), retaining
all admitted water poses in a common menu for both metals. Report every pose
and inner→outer paired transfer contrast for each source/orientation:
[(E_Dy,in−E_Dy,out)−(E_La,in−E_La,out)]. Outer-shell water is NOT a bulk chemical
potential reference. Never select raw minima across dry and wet inventories.
No hydration population, affinity or Mex–LBT3 classifier follows from this panel.
No automatic additional poses, torsion expansion or water-count sweep.

## Runnable interface

`python LBT3_WATER_PREPARE.py --config CONFIG.json --output FRESH_DIRECTORY`

Config requires `repository` plus pins `implementation`, `plan`,
`preparation_result`, `water_template`, and `dependencies` containing
`hydration_basin_coordinates.py` and `affordable_common.py`. Optional `pins`
are verified too. Require a Slurm allocation for all geometric calculations.
Runtime dependencies are Python, NumPy and SciPy; no OpenMM or energy model.
