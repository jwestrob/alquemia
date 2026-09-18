# Test native CHELPG sampling stability on the same saved densities

The verified GGR source-self shift is dominated by the projected charge vector,
not the change in source distances or effective Born radii. Determine whether
the charges and resulting source-self term are numerically stable before
changing the physical solvent model. This is a diagnostic; no new biological
score, calibration or threshold is produced. Jacob's autonomous goal authorizes it.

## Native route and qualification

Installed ORCA6.1.1 standalone `orca_chelpg` accepts GBW and optional density,
without exposed sampling flags. The native main program documents `!CHELPG`
and `%chelpg GRID/RMAX/VDWRADII/DIPOLE`. Its restart manual documents
`MORead`, `%moinp`, and `CalcGuessEnergy NoIter` for properties/energy of fixed
orbitals. Do not modify GBW binary metadata or silently request a new SCF.

First execute two property replays, Ca/La on consumed GGR2FVY, using exactly
the saved responsive orbitals, coordinates, method/basis/ECP, generating field,
charge and multiplicity. Replace EnGrad with SP; explicitly request MORead,
CalcGuessEnergy, NoIter and default CHELPG GRID0.3A/RMAX2.8A/COSMO/no dipole
constraint. Preserve every imported and resulting density artifact.

Qualification requires normal native termination, explicit fixed-orbital mode,
no SCF optimization or gradient, unchanged basis/ECP/electron inventory, default
charges matching the archived utility within2e-6e, and the replayed energy
matching the archived same-Hamiltonian energy within1e-7Hartree. Query the new
saved density at the original exterior ESP validation positions and compare
native potentials within1e-8au. A failed identity gate prevents the sampling
extension. Actual differences must be diagnosed, not called sampling effects.
The existing generic runner expects SCF convergence; retain that raw receipt
flag and use a separate explicit NoIter property-status contract, never present
a nonconverged SCF as converged. All program starts count toward measured cost.

## Conditional sampling extension, frozen before outputs

After qualification, use all five original GGR/alpha geometries (GGR1GLGextended,
2FW0,2FVY, alpha1F6S,alpha6IP9), both endpoints. Four sampling settings:

- original: GRID0.3A, RMAX2.8A;
- finer: GRID0.2A, RMAX2.8A;
- finest: GRID0.15A, RMAX2.8A;
- extent: GRID0.3A, RMAX3.5A.

All use COSMO radii and no dipole constraint. Forty native property starts in
total, including the initial two; no added SCF, new density optimization,
geometry change, solver, MACE or force call. Recheck the unchanged density
potential for every result at the same validation positions. The forty default/
variant queries are actual utility calls, not fabricated identity evidence.
Do not reuse a default result unless source/state and executable pins match.

Reuse the exact source-to-physical projection and native GK matrix for each
geometry. Report charge sums/changes, native exterior-potential quality, and
the source-self endpoint, Ca-minus-La and between-structure changes for every
setting. Preserve raw charges and all invalid cases. No per-case selection,
refitting regularization, added charge scheme or changed biological decision.

Freeze a0.5kcal/mol sensitivity screen for each paired source-self contrast
and each of the three GGR structural differences: finer-to-finest and original-
to-extent changes must each be within0.5kcal. This is below the prior1kcal
source-representation energy screen and far below the18kcal discrepancy; it
is not a classification threshold or proof of asymptotic convergence. Also
report original-to-finer changes and all endpoint effects without hiding them.
Potential/charge quality gates remain those of the parent observation protocol.

Use the existing ORCA runner and MPI/thread fixes, with32CPU for the initial
two16-rank replays and64CPU/up to four concurrent16-rank replays if extended.
Density utilities use one CPU each in a separate appropriately sized allocation.
No project CPU/time budget. This is a finite declared diagnostic inventory;
no automatic production refit or rescore. The original10SCF stage cost426wall
seconds on64CPUs; NoIter cost must be measured and is not assumed free.

Sources: [ORCA population analysis](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html),
[ORCA fixed-orbital replay](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html).
Manuals and standalone help are preserved under `chelpg_sampling_*` workspaces.
