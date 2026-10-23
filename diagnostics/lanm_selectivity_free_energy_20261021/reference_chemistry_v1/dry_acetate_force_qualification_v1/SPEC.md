# Matched dry mono-acetate: sampled reference derivative qualification

Declared before displaced-coordinate preparation or molecular execution. Root
owns release/submission. This is a consumed-reference numerical diagnostic,
independent of the running water fit. It opens no reserved validation outcome.

## Question, existing coverage and exact centers

Do the saved uncorrected PBE0-D4 forces for the exact dry mono-acetate La/Dy
centers agree with energy derivatives for rigid radial translation and ligand
rotation? The consumed donor-fit force-direction audit found large dry donor
discrepancies, but this does not establish which model/reference term is wrong.

The searched project source/declarations contain water-dimer radial qualification
and central **hydrated** acetamide translation qualification. Neither covers any
of the six matched-dry donor centers or ligand rotations. Matched-dry CP/CC
studies compare energies; planned model testgrad checks qualify a different
Hamiltonian. Native analytic-gradient admission alone is not a derivative test.
All six dry authorities and this limited coverage finding are retained in INPUTS.json.

Use training99ea27f8 and its exact `matched_dry_acetate_mono_central_La/Dy`
rows. Central reference RESULT hashes are4301f2f0 and b4f813d3; parent native
panel304da68b. Eight real atoms, charge+2, effective singlet,40 explicit
electrons,189 orbital functions; physical La1/Dy6 and La46/Dy55 cores remain
separate. Metal index0 and intact ligand indices1–7; coordinating O3, other O4
(zero-based). Reuse central energy and full gradient. Do not recalculate centers.

## Exactly sixteen new single-point energies

Two metals × two geometry-defined directions × two steps × two signs.
Set q=−0.010,−0.005,+0.005,+0.010Å, keeping the metal fixed:

- Radial: translate every ligand atom by q u, using the stored outward unit
  direction u. It equals normalized metal-to-coordinating-O3.
- Rotation: let c be the arithmetic centroid of all seven ligand atoms and
  n=normalize[u×(R_O4−R_O3)]. Freeze c and n from the original geometry;
  rotate the complete ligand about c with right-handed Rodrigues angle
  theta=q/L, **L=1Å**. This is the existing torque/1Å coordinate convention,
  not an outcome-selected torque direction. The two angles are0.010/0.005rad.
  A degenerate axis is a preparation failure; no fallback/search is authorized.

Both metals use identical coordinates and motions. Serialize12-decimal XYZ;
verify source labels/order, fixed metal, exact prescribed transformation within
1e−12Å, and rigid internal distances within1e−9Å. These are inherited geometric
precision checks. Preserve all16 IDs and every failed/missing slot.

## Recipe and reused implementation

Reuse the original executed center deck and assets via original REFERENCE.py
eec0e60d, the hydrated-amide thin-adapter pattern7de9c17a, and Water48 strict SP
parser2871d929. Preparation must replay both original raw results and exactly
reconstruct their executed input at recorded resource settings. For each new
SP, remove only `EnGrad`; geometry changes as above and PAL/MaxCore follow the
actual allocation. Keep PBE0 D4 def2-TZVP def2/J RIJCOSX, NoAutostart, DefGrid3,
VeryTightSCF, PModel, MaxIter300 and exact external orbital/ECP/Aux-J bases.
No restart orbitals, CP, CC, solver/grid changes, new gradients or model calls.

Reuse native state/coordinates/basis/electron/ECP/termination/SCF-stop and
electronic+D4+final-energy closure admission. Preserve printed residual flags,
active stopping mode and final COSX/integration metadata. The center is EnGrad
and displacements are SP by design: stage comparison is reported, not assumed
to prove a common numerical surface.

## Signs, units and prospective resolution

At each h, d_FD=[E(+h)−E(−h)]/(2h). Report electronic, D4 and total derivatives
separately, using the original627.509474kcal/mol/Eh conversion. Saved total
analytic derivatives are:

- radial: −sum_i F_i·u over ligand indices1–7;
- rotation: −n·sum_i[(R_i−c)×F_i]/L, with L=1Å.

Thus all reported coordinate derivatives have units kcal/mol/Å. The opposite
sign is generalized force. Report La, Dy and **Dy−La** independently for both
directions. Translation net-force closure remains visible; do not replace the
ligand derivative with a metal projection.

Use the existing **0.1kcal/mol/Å** diagnostic resolution, pinned to Water48 SPEC
f313f90a and hydrated-amide SPEC f5de75cf. For each of the six direction/comparison
channels require total FD(h)−analytic within0.1 at **both** steps and the two
total FD estimates within0.1:18 primary consistency checks. Retain electronic
and D4 two-step differences/0.1-resolution flags separately. Use the0.005 result
as the prespecified summary only when the corresponding two-step check passes.
No fitted tolerance, extrapolation, third step or removal of an unfavorable row.

Agreement qualifies only these two sampled directions of the two exact centers
at this resolution. A resolved FD/analytic discrepancy is sampled numerical
inconsistency; unresolved step dependence can reflect curvature or numerical
noise. Neither identifies SCF failure or a unique donor-physics defect. No
physical-state/basis/method accuracy or other dry donor direction is certified.

## Finite staging, failures and costs

First release exactly the four **La radial** SPs. Preserve their results in the
16-cell inventory. Remaining12 may run only after actual raw/input admission,
the three La-radial primary checks and a pinned root scout review pass. Reuse
those four original outputs; no duplicate scout or central call. Root separately
binds exact stage source, panel, output directories and finite result inventory.

Each singleton uses4MPI×1thread, standard_expand, scheduler-default memory/time
and MaxCore=floor(0.75×actual granted MiB/4). No automatic retry or requeue.
Record wall time and allocated CPU separately; do not infer accuracy from cost.
Saved preparation/collection use one local CPU and make no native calls. A failed
or absent cell remains unavailable; complete native counts cannot override a
failed numerical check. All16 are required for the full two-direction result.
