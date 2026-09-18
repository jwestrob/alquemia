# Prepare a charge-closed density/GK hybrid boundary

Declared while the independent permanent-multipole density utilities are running.
This step prepares physical states and native parameters only. No scored hybrid,
new energy, response solve, DFT, MACE or geometry operation is part of this step.
Active-goal autonomy applies. The production baseline stays unchanged.

## Fixed real systems and source definitions

Use GGR extended/connected and alpha1F6S/6IP9 with the existing normalized
physical coordinates, peptide-amide cores, cap projection, protonation,
disulfides, waters, assembly and evidence strata. Restore the actual physical
metal to the native framework inventory:4698/1932/1898atoms. Its position is
the recorded source-metal position, not an inferred or optimized site.

Keep the same QM source support including all physical cap-projection anchors.
Use the actual projected CHELPG distributions from normalized_charge_v1 only
as frozen reaction-field sources. Do not replace them with the failed fitted
charge/dipole models. Preserve printed fit charge residuals; no renormalization.
Ca/La source formal charges are−1/0 with the same−3ligand formal-charge ledger.

## Explicit environment-charge boundary

Read the full actual AMOEBA2018 framework monopoles/multipoles. For every
residue intersecting Q support:

1. Derive its formal charge from the complete AMOEBA residue sum, requiring
   agreement with an integer within1e−8e and with the existing chemical ledger.
2. The exterior target charge is that formal charge minus the recorded
   selected-fragment formal charge. Extra cap-anchor residues have selected
   formal charge zero. Completely selected neutral waters retain no exterior.
3. Remove all FF permanent moments on Q-support atoms. Starting with the actual
   remaining AMOEBA monopoles, distribute the residue's required residual equally
   over the distinct exterior atoms directly bonded to Q-support atoms in that
   same residue. Use the source graph, not sequential residue numbering.
4. A nonzero residual without eligible recipients is unsupported. A fully
   selected residue must have target zero. Do not borrow remote recipients,
   omit a residue, change protonation, or globally neutralize the protein.

This is a new explicit AMOEBA boundary rule. Reuse existing graph/fragment
identities, but do not copy ff19SB partial charges or numerical increments.
Retain every original/modified charge, affected residue, recipient and bond.
Exterior dipoles/quadrupoles remain the actual native tensors. All physical
axis atoms and the complete covalent graph stay present even where moments are
removed. Apply source-moment replacement once, after native parameter matching.
Freeze induced variables on every Q-support atom independently of the source
moments. Environmental polarizabilities/damping remain native AMOEBA2018.

Pinned `kpolar` removes a zero-moment, zero-polarizability source metal from
its active multipole index list. After source replacement the adapter must
explicitly restore the full physical `ipole/pollist` mapping, including that
source, and its recorded multipole order. This restores indexing for the
declared source moments; it must not add a fabricated initialization charge.
Export the actual mapping. Native polarization groups are already constructed
over all physical atoms and must remain unchanged. Zero-source states retain
the same mapping with zero source moments, rather than changing the cavity or
losing an atom from later reaction-field sums.

Expect exterior formal charge−3for GGR and−4for alpha; require actual sum within
1e−9e. Check per-residue closure rather than accepting a coincidental whole sum.
The paired permanent environments, coordinates, maps and permitted induced
variables must be identical between Ca/La, not merely have equal net charge.

## A common physical source-metal cavity, explicitly not La force-field fitting

Use the pinned AMOEBA2018 Ca type358/class99 cavity specification for **both**
source endpoints. Declare separate local QM-source atom types with actual
atomic numbers/masses20/40.078 and57/138.90547. Give both the same class99
vdW/cavity support and SOLUTE tuple3.6670,3.6670,3.6497,0.1350. These are the
actual native parameter-file values, not fitted to a result. The source-metal
permanent moments are initially zero and replaced only with the QM-projected
monopole; classical source-metal polarizability is zero and response disabled.
This is **not** an AMOEBA La ion parameterization, a formal+3point ion, or use
of La09/La22 response parameters. Native vdW/bonded energies remain disabled.

Pinned `ksolv.f` predicts source-metal GK radius1.82485Å, descreen radius1.795Å,
overlap factor0.72 and neck factor1for this unbonded heavy source. The SOLUTE
neck parameter is superseded by native nheavy=0 behavior; record the actual
values rather than asserting0.135is the final neck factor. Keep the prior
descreen offset0.30Å, dielectric offset0.09Å, internal dielectric1, native
bulk78.3, GKC2.455, Grycuk/neck/tanh settings. Report all cavity/dispersion
parameters and verify they match across metal labels; do not silently alter
native constants. Nonpolar terms will cancel only after matching is proven.
Retain the already qualified response settings1e−7Debye and100iterations in
these native keys, although this preparation executes no induced solve.

Class99 is used for the common cavity/dispersion support, not to add a Ca
short-range metal interaction into the MACE hybrid. This common source cavity
is a declared approximation. Its radius sensitivity belongs in the future
hybrid validation, not an adaptive change after inspecting classifications.

## Native preparation inventory and acceptance

Four representations × four parameter states =16native initializations:
Ca source, La source, zero-source/Ca identity, zero-source/La identity. Each
exports parameters and actual local/global moments only; no Born calculation,
energy(), induction, force or atom motion. The two zero-source states establish
whether one environment-only reference may later be shared for a pair. They
may share a reference only after the complete actual cavity, permanent and
response parameters are verified identical apart from inert element/mass labels.

Use a separate versioned preparation and isolated frontend; preserve existing
native source/library and all experiment outputs. Source/context hashes and
full source-to-native atom order are recorded. Fail explicitly for missing
cofactors, unsupported residue charges, absent recipients, source/axis mismatch,
charge/parity inconsistency or changed water/assembly/protonation state.

Checks: exact paired physical coordinates; source-moment replacement and zero
source-induced-variable masks; moment local/global reconstruction≤1e−10native
units; cavity equality and common declared metal parameters; same source graph;
exterior closure1e−9e; projected QM sum within1e−5e of formal charge, with the
actual residual reported and preserved; full endpoint charge equals exterior
plus actual projected QM sum. Do not use rounded charge totals in energies.

Expected scale: existing three native moment initializations cost~2.7summed
wall-seconds;16initializations are short preparation work, not a new production
cost regime. Run each with one native thread and record actual costs/failures.
No project resource/time budget. Prepare the subsequent finite energy/field/
response manifest only after the boundary and density-coupling prerequisites
qualify. That later model retains the existing partition and evidence rules.
