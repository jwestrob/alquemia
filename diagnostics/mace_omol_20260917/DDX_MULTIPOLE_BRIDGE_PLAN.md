# Qualify the full frozen multipole representation before a transfer score

The native frozen-response functional passes on four GGR endpoints. Its
conductor counterpart needs both the surface potential and source integral
for the same permanent/induced distribution. Qualify that representation
independently of the running discretization test; do not compute a new score.

Use the archived GGR 2FW0/2FVY, primary/rigid, Ca/La static multipoles and tight
supplied-density response dipoles from collection_job_1201162. For each of
these eight real states, form exactly the two algebraic distributions already
specified in the energy derivation:

    average: P + (mu_d + mu_p)/2
    difference: (mu_d - mu_p)/2.

The source atoms have zero induced dipoles in both. These are algebraic
distributions, not new molecular microstates. Preserve every physical cavity
sphere, intrinsic radius, atom order, source charge and permanent environmental
moment. No charge renormalization or force-field charge on source atoms.

## Native interface and matching representation

The stock ddX 0.9.0 module uses real spherical multipoles ordered by l=0,1,2
and m=-l...l. Verify normalization/signs against its actual harmonic and
multipole source routines. Convert the native Tinker global moments to atomic
units once: dipoles divide by BOHR_TO_A, quadrupoles by BOHR_TO_A squared.
For q, Cartesian dipole d and native effective Cartesian quadrupole Q:

    M00 = q/sqrt(4*pi)
    M1 = sqrt(3/(4*pi)) * [d_y,d_z,d_x]
    M2 = sqrt(15/(4*pi)) * [2Qxy,2Qyz,
                  (2Qzz-Qxx-Qyy)/sqrt(3),2Qxz,Qxx-Qyy].

Native Tinker couples its effective Q directly to the potential Hessian;
do not apply its original parameter-file factor of 1/3 again. The actual
rounded moments have a small nonzero trace (up to about 9.4e-7 e*A^2 observed
in the first real state). Record it. For a point multipole outside its center,
the trace cancels because the Coulomb potential is harmonic. Do not alter the
archived GK moments or pretend the stored trace is exactly zero.

Use four fresh native Model setups, one per physical geometry, lmax=6,
Lebedev=194, conductor/epsilon78.3/eta.1/shift0, 64 threads, no forces and
enable_fmm=False for these source evaluations. This explicitly obtains dense
native source potentials; no forward or adjoint continuum equation is solved.
Evaluate sixteen native multipole potentials and sixteen native psi arrays.
Retain the full arrays and all source/build/runtime hashes.

## Frozen checks

- Native dense potential versus an independent Cartesian Coulomb, dipole and
  quadrupole evaluation at 256 evenly distributed actual cavity points:
  maximum difference <=1e-10 atomic units. This samples the full distribution
  at declared real points; do not describe it as exhaustive spatial coverage.
- Native psi versus `4*pi*M_lm/[(2l+1)*radius^l]`, including zero higher modes:
  maximum error <=1e-12 in native units.
- For each distribution, evaluate the primary probe points transformed by the
  archived rigid rotation/translation using the actual rigid-state moments.
  Potential change <=1e-8 au; record input moment/dipole transformation errors.
- All arrays finite; native input parameters and physical source coordinates,
  radii, identifiers and source frozen-dipole inventory match their archives.
- Record actual charge sums for each distribution. The difference distribution
  has identically zero charges; no invented total-protein neutralization.

Four native setups/16 source-property calls/16 source integrals, zero continuum
solves, response iterations, quantum calls, MACE calls, force calls or scores.
One 64-CPU/128-GiB allocation. Measure costs; no automatic retry or project
compute/time cap. Passing this interface test does not qualify grid convergence,
self-consistent solvent response or biological predictive value.
