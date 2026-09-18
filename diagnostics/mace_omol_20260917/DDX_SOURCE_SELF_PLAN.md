# Resolve the protein dielectric boundary for the fixed source-self term

## Question and scope

Does the verified GK source-self discrepancy between GGR 2FW0 and 2FVY persist
with a resolved PCM boundary, on the same nuclear geometry and fixed default
projected charges? Sampling alone left a 14.97–18.13 kcal discrepancy. This
experiment examines the solvent approximation separately. It is development
on two consumed structures of one protein, not a biological validation panel.
The autonomous discriminator goal authorizes execution.

Use the actual primary and archived rigid-transform static outputs from
`trial_gk_expansion_native_v1/collection_job_1201162.json`. Preserve the complete
physical atom inventory, coordinates, native GK intrinsic radii (`rsolv`,
exported as `radius_A`), source support and default projected CHELPG charges.
Both metals use the same physical cavity, including the existing common
1.82485 A source-metal radius. All exterior atoms remain as cavity spheres;
their source charges are explicitly zero for this source-only diagnostic.
This is not the environmental identity limit and does not model the full
protein electrostatic energy. No water, cap, protonation or donor changes.

## Model and energy

Use the isolated, pinned ddX 0.9.0 PCM solver. Interior relative permittivity 1,
exterior 78.3, salt zero; union of the specified atom-centered spheres, smoothing
eta 0.1 and shift 0. No nonpolar term, induced polarization, direct protein
coupling, CPCM quantum term, MACE inference, gradient or structural response.
The selected GK radii become an explicit PCM cavity specification. GK's neck,
descreen and tanh approximation is not asserted equivalent to this boundary.
Agreement with GK is not an acceptance criterion or a target for tuning radii.

Both surface potential phi and density integral psi come from the same
fixed atom-centered source monopoles through ddX's multipole interface.
`E_self = 0.5 * dot(psi, x)` Hartree, where x solves the PCM equations.
Use the library energy and independently audit this contraction. Convert once
with the existing 627.509474 kcal/mol/Hartree constant. Report Ca-minus-La
and 2FW0-minus-2FVY component contrasts, plus their differences from native GK.
No full hybrid score, aquo reference, threshold or biological classification.
This remains a projected-charge model; it is not exact quantum-density PCM.

## Frozen numerical inventory and checks

Matrix-free operator; FMM enabled with multipole/local orders 12/12 at all
resolutions. DIIS history 20, solver tolerance 1e-10, maximum 300 iterations;
nonconvergence is an explicit numerical failure. Disable force machinery.
Use 64 OpenMP threads per solver and one solver at a time, single-thread BLAS.

Eight model setups and 20 forward solves:

- Both structures, Ca/La, at lmax/grid 6/194, 9/302, 12/590: 12 solves.
- Both archived rigid structures, Ca/La, at the primary 9/302 setting: 4 solves.
- One zero-source solve per original physical cavity at primary settings: 2.
- Repeat both 2FVY primary endpoints from fresh solver states: 2.

Every setup preserves paired radii/coordinates and source charge inventories.
Audit native source potential against direct Coulomb sums at every exposed
cavity point: maximum error <=1e-8 au. Verify monopole psi normalization and
zero higher moments within 1e-12. No charge renormalization.

Freeze acceptance before outputs: primary-to-refined changes <=0.1 kcal for
each endpoint, each Ca-minus-La contrast and the between-structure contrast;
rigid-transform changes <=0.05 kcal for the same quantities. Report coarse-to-
primary changes separately. Repeat and zero-source energy error <=1e-8 kcal;
library versus direct energy contraction <=1e-8 kcal. Source cross-contraction
reciprocity <=0.05 kcal for each paired model. Passive reaction energies must
be nonpositive within 1e-8 kcal. Retain every failed check and unscorable row.
These tolerances are below the existing 0.5 kcal representation screen and
well below the 15–18 kcal signal being examined; they are not fitted to a label.

Request one 64-CPU, 128-GiB allocation; record setup, potential, solve, total
wall/CPU, memory, iterations and all failures. No project time/compute stopping
budget and no automatic scientific retries. The eight model setups are serial;
free each before the next. This one-time convergence inventory is distinct
from production cost, which remains unestablished. Preserve all input mappings,
software/source pins, actual solver coefficients, potentials and receipts.
