# Mex single-proposal physical response pilot

Declared after the complete12709 directional derivative check, before computing
new proposals or energies. Preserve static12685 and all16 derivative cells.
This is one bounded response step from the same repaired Mex8FNS origin for each
metal, followed by paired cross-scoring; it is not an optimization campaign.

## Unchanged physical model

Use the admitted2049-atom physical model,209 local atoms including four caps,
1840 shifted external charges, fixed threeNd(III) spectators and charge+2 global
unshifted FF/TIP3P/formal3 vector. Score exactly
`MACE_vac(Q) + GFN2_emb(Q;MM) - GFN2_vac(Q) + MM_ledger + ALPB(qFF,R)`.
Local physical states are charge0/La multiplicity1 or Dy6; GFN2 retains its
published effectiveUHF0 convention. Preserve model, parameters, runtime,
cap rules, classical exclusions, solvent and all native admission gates.

## Frozen source graph and seven coordinates

Read the exact source export bonds and physical source IDs already bound by12709.
The four axes, in ancestor-first order, are:

1. `A/83//N` to `A/83//CA`: Leu83 backbone hinge.
2. `A/84//CA` to `A/84//CB`: complete Asp84 side chain.
3. `A/86//CA` to `A/86//CB`: complete Asp86 side chain.
4. `A/95//CA` to `A/95//CB`: complete Glu95 side chain.

For each axis remove that one native bond and take the complete connected
component containing the second atom. Require a graph bridge, membership only
in source protein atoms, correct residue identities, and nested or disjoint
ancestor-before-descendant membership. No truncated segment, individual donor
motion, synthetic cap coordinate or alternate axis may substitute after failure.
The last three coordinates are Cartesian translations of the EF3 target alone.
Waters and all spectator ions remain fixed. Save exact source graph, axis/moving
IDs, physical order, origin Jacobian and all metric values before energy calls.
Require each acidic side chain to retain its full carboxylate group; protect all
source carboxylate/amide groups, aromatic rings and actual peptide groups (with
flanking CA and N-bound H) by their internal pair distances. After excluding an
axis's two endpoints, no protected group may be split by its moving component.
Original peptide and chirality guards cover the complete source protein.

## One own-gradient proposal per metal

Use the archived total physical gradient of each own metal. Let `J` be the
6147-by7 origin physical-coordinate Jacobian. Normalize each nonzero column by
its Euclidean norm, forming `U`. Retain eigenvalues of `U.T@U` greater than
`1e-10` times its largest eigenvalue and use that Gram pseudoinverse. Form the
negative own-gradient coefficients, then convert back to radians/angstroms.
Do not use the La/Dy difference or experimental labels to choose a direction.

Apply exactly one common scalar to each proposed direction. The scalar is the
minimum of (a)0.02Å divided by the largest per-atom sum of absolute torsional
arc-length and translation bounds and (b)0.5degrees divided by the largest
absolute torsion coefficient. This is the same conservative source-kinematics
rule as the earlier PREPARE.py. Compose hinge before descendant rotations using
current axes and Rodrigues rotations, then target translation. Require negative
own linearized total work. No backtracking, alternate basis, retries or step-size
adjustment after inspecting geometry or energy.

Require all-physical displacement ≤0.02Å, each angle ≤0.5degrees, fixed-atom
identity, source bond change ≤1e-10Å, angle-cosine change ≤1e-10 and protected-group
pair-distance change ≤1e-9Å. Reuse original source chirality/peptide-basin checks.
Check severe contacts for every moving atom, including the current moving target,
against all2049 atoms relative to the immutable original source, excluding native
1–2/1–3 neighbors: reject distances below `min(source distance,1.0Å heavy or
0.55Å involving H)−1e-7Å`. Existing source short contacts remain recorded.
An unavailable proposal remains an unavailable declared slot; no replacement.

## Coherent coordinates and finite matrix

For each admitted proposal update every real QM atom, MM field position,
zero-field boundary anchor, native/export physical record and boundary-coordinate
record from the same physical array. Keep charge-shift membership and every
charge/hardness fixed. Reconstruct all four caps and both anchor Jacobians from
current physical anchors; unchanged anchors preserve their already validated
coordinates. No cap enters the physical solvent cavity. Use round-trip17digit
coordinates and verify source-ID/element/charge inventories and both metals'
identical physical coordinates before the first new molecular call.

Evaluate each admitted proposal with both metals. The maximum is four new local
MACE forwards, eight local GFN2 calls (embedded and isolated), four complete MM
configurations/sixteen group queries, and four full fixed-charge ALPB calls.
Reuse the two exact origins. Use the existing pinned CPU Python/Reference OpenMM
subprocess route and warm native MACE; at most four one-thread GFN2 workers.
Native MACE may return its existing analytic forces in the same forward; these
are retained as receipts and do not trigger further steps. Request MM forces in
the same sixteen group queries, then retain complete candidate total/component
gradients and current-axis loads scaled by the original Jacobian column norms.
This uses already returned native gradients, adds no calls, and supplies a
residual for any separately declared later response decision. No new DFT, whole-
protein SCC, minimization or parameter fit. Root owns Slurm submissions and costs.

## Report and stopping rule

Preserve all four proposed state/cell slots, exact provenance and actual attempts.
Report both metals' predicted own linear work, actual own component/total work,
cross work, common strain and Dy−La work for each geometry. Form the shared menu
`{origin, La proposal, Dy proposal}` only when both proposals and all four new
cells are available. Both metals minimize over that identical menu; report
selected geometries and the change in their conditional Dy−La contrast from the
origin. Difference native component energies before conversion/summation to
avoid subtracting unrelated large total-energy references.

Own descent means actual own total work ≤0 (with no fitted allowance). Geometry,
execution or own-descent failure closes this finite pilot and prevents expansion;
valid cells and the origin remain separately available. Do not replace a failed
proposal with the origin and claim a complete common3 test. A successful tiny
feasible response is not a stationary minimum, a full-protein affinity label or a
Mex-versus-LBT comparison; EF3 remains conditional on three fixed Nd spectators.
