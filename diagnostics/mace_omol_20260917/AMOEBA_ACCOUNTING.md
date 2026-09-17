# Polarizable-environment accounting and unresolved capabilities

This is a design/accounting result, **not an implemented scoring protocol**.
No AMOEBA energy/force has yet been evaluated on these cases.

## Maintained solver behavior verified

Installed OpenMM 8.5.1 Python builders were inspected and pinned. The matching
upstream 8.5.1 reference implementation is preserved in
`workspaces/mace_omol_20260917/amoeba_sources_v1/`. The installed multipole
header is independently recorded in `verification.json`.

- `AmoebaReferenceKernels.cpp`, multipole setup/execution: when GK is present,
  the multipole kernel constructs the GK-coupled implementation and supplies
  its dielectric, radii, charge, cavity and descreening parameters.
- `AmoebaReferenceMultipoleForce.cpp`, GK `calculateInducedDipoles`: vacuum
  and solvent fields enter coupled induced-dipole solutions. The solvent
  changes the polarization; this is not a fixed-dipole GK add-on.
- GK `calculateElectrostatic` includes vacuum electrostatics, all GK pairs
  including diagonal reaction/self terms, the enabled cavity term, and the
  vacuum-to-solvent response correction. The GK force's separate reference
  `execute` returns zero because the multipole kernel handles its terms.
- Consequently, a force-group query for GK alone does not isolate a usable
  solvent energy. Matched separate Hamiltonian evaluations are needed to
  decompose permanent, induction and solvent effects. Reference-platform
  source inspection does not replace future CUDA integration checks.
- Permanent, direct-polarization and mutual-polarization covalent scales are
  distinct. The framework stores all eight maps. The AMOEBA XML advertises
  permanent 1–4/1–5 scales 0.4/0.8 and separate polarization-group rules.
  They cannot be replaced by unscaled Coulomb coupling without a new model.

Primary source links:
[multipole implementation](https://github.com/openmm/openmm/blob/8.5.1/plugins/amoeba/platforms/reference/src/SimTKReference/AmoebaReferenceMultipoleForce.cpp),
[reference kernels](https://github.com/openmm/openmm/blob/8.5.1/plugins/amoeba/platforms/reference/src/AmoebaReferenceKernels.cpp),
[multipole API](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.AmoebaMultipoleForce.html).

## A concrete subtraction target, conditional on a valid source model

Let Q be the QM distribution projected onto physical source atoms; E the
remaining permanent protein multipoles; A the permitted environmental induced
dipoles; and B the **same full physical cavity** in every subtraction below.
Core/source atoms have no independent induced variables. Local-axis anchors
remain physical even when their own permanent multipoles are replaced.

Define F_B(Q,E; A_E) as the maintained solver's converged permanent plus
induction plus GK energy, with all physical positions, radii, damping and
covalent maps explicitly fixed. Then the environmental transfer target is

```
DeltaU(Q) = F_B(Q,E; A_E) - F_B(0,E; A_E) - C_vac(Q,Q)
E_trial(M) = E_DFT,vac(M) + DeltaU(Q_M)
             + T_MACE,short(full,M) - T_MACE,short(core,M)
```

`C_vac(Q,Q)` must be computed with the **same source-source permanent scaling**
used by F, with environment permanent multipoles zero, induction disabled and
GK absent. It removes the classical core electrostatics already represented
by DFT. It is not the isolated quantum energy or an unscaled arbitrary pair sum.

`F_B(0,E; A_E)` retains all source atoms as cavity/axis sites and removes their
permanent moments, keeping the same environmental response model. It cancels
environment-only permanent/induction/solvent terms. It is identical between
Ca/La endpoints only if that fact is verified. It need not be identical across
different numerical QM partitions and must never be reused on identity alone.

The difference retains direct Q–E interactions with declared boundary scales,
environmental induction caused by Q (including its coupling to solvent), and
the change in reaction field, including Q self terms. Geometric cavity terms
cancel only for the identical B/parameters; exclude bonded/vdW/WCA terms from
these electrostatic F evaluations. No CPCM term is added: the proposed parent
is the already computed vacuum DFT state. A CPCM parent needs a different,
matched transfer definition. An embedded-responsive DFT parent also needs
its external-field expectation handled explicitly; it cannot be substituted
into this expression while retaining all of DeltaU.

The identity limit removes E and its polarizabilities **and** makes B the
vacuum reference, leaving C_vac which cancels. Merely zeroing protein charges
does not satisfy that limit.

## What remains unresolved before this can produce a score

1. **Source damping.** OpenMM's per-atom minimum rule does not implement the
   La paper's POLPAIR model. Zero polarizability also sets the default damping
   length to zero through its builder, removing damping in the reference
   field expression. Thus “set core alpha to zero” is insufficient. A supported
   frozen-source scheme must preserve justified finite damping independently
   of the allowed induced variables. Tinker's `POLARIZABLE`/POLPAIR facilities
   are the next capability target, not a proven solution.
2. **Boundary multipoles.** Replace FF permanent moments on *all* projected-Q
   support sites, including cap anchors, exactly once. Keep source identity,
   local-axis anchors and the full graph. The old ff19SB monopole redistribution
   cannot silently become an AMOEBA multipole rule. A new explicit residue
   formal-charge ledger and boundary redistribution rule must be declared and
   tested. There is no approved numerical hybrid boundary in this stage.
3. **Coverage/units.** The three protein/water frameworks pass. La, its damping,
   GK radius and compatible Ca treatment are not parameterized by that pass.
   PQQ/nonstandard-cofactor coverage has not been tested. Published La09/La22
   parameters cannot simply be relabeled AMOEBA2018.
4. **Field accuracy.** Earlier potential/coupling checks do not verify the
   electric fields driving induction. Compare the actual native QM field with
   the projected representation at environmental sites before interpreting
   induction energies. No such field test has run here.
5. **MACE double counting.** The qualified short readout is jointly trained and
   is not a uniquely separable non-electrostatic energy. Adding the above
   correction is an explicit hybrid approximation; success must survive
   physical/partition and predictive checks, not just algebra.
6. **Acceptance/cost.** Retain the preexisting 2 kcal/mol partition and 0.01
   numerical energy gates; separately declare induction convergence and
   field-quality tolerances before their results. Numerical score, absolute
   reference, classification, relaxation, entropy and combined gradient remain
   unavailable. Preparation receipts do not measure solver production cost.

Parameter arrays in `mapping.json` use OpenMM MD units. Multipole order is
charge (e), dipole (e nm), quadrupole (e nm²), axis enum, Z/X/Y anchor indices,
Thole coefficient, damping factor (nm^(1/2)), polarizability (nm³). GK order is
charge (e), radius (nm), scale, descreen radius (nm), neck factor. The complete
serialized systems retain original typed quantities and force parameters.
