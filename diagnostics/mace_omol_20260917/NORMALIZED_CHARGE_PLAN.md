# Qualify saved normalized-core charges for a full-boundary solvent model

Declared after MATCHED_H_REPORT.md, before charge utility outputs. Active-goal
blanket authorization applies. This is a prerequisite for the solvent extension,
not a new affinity score and not a modification of the failed vacuum result.

## Question, inputs and method

Can one uniform, inexpensive QM charge representation reproduce the potential
of the eight normalized vacuum endpoints, and can their artificial cap charges
be mapped to actual protein atoms without a large additional error?

Use only completed job1200950, manifest `matched_H_quantum_v1/manifest.json`,
and its exact source mapping in `matched_H_prepared_v1/preparation.json`.
Four consumed representations: GGR extended/connected, alpha1F6S/6IP9; Ca/La
pairs. No new SCF, geometry, donor, protonation, water or assembly change.
The native GBW, densities and densitiesinfo files survive. Copy them into
isolated task directories, pin their hashes, and preserve originals.

Eight `orca_chelpg saved.gbw` calls with ORCA6.1.1 native defaults:0.3A grid,
2.8A extent,COSMO exclusion radii,no dipole constraint. Verify output echoes,
atom order, total charge and native ECP convention. No charge-scheme comparison
or selection by protein. Eight `orca_vpot` calls on the same saved densities
provide independent potentials at the established exterior probe set:
32 Fibonacci directions per atom,1.4/1.8 radius shells, excluding all atom
interiors at1.4 radii, using the existing affordable_solver.real_esp_points.
Common Ca/La sampling radius1.8A; these are potential probes, not a solvent cavity.
Ca/La points must match exactly. Charges are parsed at native printed precision;
no renormalization or invented unprinted digits.

The native documentation describes total-charge-constrained ESP fitting and
allows standalone evaluation from GBW. We checked the6.1 manual and existing
successful utility implementation. This repeats that established scheme on
new compatible wavefunctions; CHELPG itself is not a new development.

## Cap mapping diagnostic

Build a fixed linear projection from the source graph. A source atom charge
maps identically to its real atom. For a sigma link H at
x_cap=(1-lambda)x_retained+lambda*x_omitted, lambda=cap_length/bond_length,
put (1-lambda)q_cap on the retained atom and lambda*q_cap on the omitted atom.
The real coordinates are unchanged. This preserves the represented monopole
and dipole exactly, but changes higher moments; quantify that error against
both original fitted and exact QM potentials. Reject missing anchors, duplicate
source identities or lambda outside(0,1). Overlapping contributions add through
one source-index matrix. No synthetic cap sphere or coincident cap charge is
introduced. No forcefield charge is added in this diagnostic.

This projection does not by itself solve classical boundary charge closure.
Its validity is assessed independently of any biological ordering. A future
full charge model must document how fixed exterior protein charges and local
covalent boundaries are combined; do not fill that missing model with zeros.

## Frozen checks and outputs

Per endpoint, require finite ordered charges, sum within5e-5e of formal charge
(native printed precision), and exterior potential RMS<=0.005au OR relative
RMS<=0.10, the same historical charge-quality gate. Apply the same rule to
projected potentials; report projection-only error separately. Projection
monopole/dipole conservation tolerance1e-9e and1e-8eA. Report matched Ca-minus-La
potential errors as well, using identical probes, with the same RMS gate.
This is a representation check; it cannot guarantee sub-kcal solvent energies.
No electrostatic energy, affinity, calibrated class or combined gradient is
reported. A later solvent-energy test still needs its own quantitative checks.

Use eight CPU workers in a shared allocation, one thread per utility,16GB total
host memory. Sixteen utility calls, no new DFT/MACE forward/solver/training.
Previous four-state CHELPG work cost296allocatedcore-seconds; larger111-atom
cores can take longer. Expected seconds to minutes, no application time/CPU
budget. Use existing utility executor/resource receipts, frozen implementation,
exclusive locks, immutable attempt directories and explicit partial recovery.

Persist parsed charges, exact/projected potentials, projection matrices, all
source/input/executable hashes, resource receipts and all failures. Tests use
these real pinned inputs; actual integration stays unrun until outputs exist.
Baseline/default remain unchanged. No automatic solver or prediction follows
a nominal charge-quality pass.

## Energy accounting for the intended extension

The candidate being investigated is H_M(vacuum)+G_full[q_M], where
H_M is the recorded matched-H DFT+MACE(full-core) expression and G_full is
reaction-field energy on a common physical full-protein boundary. It would
add no bare direct Coulomb term: learned vacuum interactions already occur
in the context term. It would add no core CPCM energy. For a fixed-boundary
linear GB model, G=-1/2(1/epsilon_in-1/epsilon_out) sum(q_i q_j/f_GB,ij),
including self and cross terms. Identical environment-only terms cancel between
Ca/La; core/self and core-environment reaction terms do not. This is a frozen
charge approximation, with no solvent-induced change of the quantum density.
The full charge model and solvent execution are not yet declared or implemented.
Do not claim that adding G repairs omitted long-range vacuum physics in MACE.

Sources: [ORCA6.1 standalone CHELPG](https://www.faccts.de/docs/orca/6.1/manual/contents/utilitiesvisualization/utilities.html#orca-chelpg),
[OpenMM GB energy](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#gbsaobcforce).
These document components, not validation of this assembled discriminator.
