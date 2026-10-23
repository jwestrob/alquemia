# Finite initialized-PME export qualification

Prepared only; linking and native execution require root binding. This tests the
export and reciprocal operator on one consumed neutral water dimer, not periodic
response stability, physical accuracy, a full box, or dynamics. All fatal
periodic-response guards stay intact. No selected candidate is used.

## Fixed fixture and settings

Use the exact six-atom `aquo9_central_La__near__two_water` input from the admitted
57-program panel. Preserve its XYZ and parameter bytes and its key prefix.
Append orthorhombic axes 44 Å, angles 90 degrees, `EWALD`, real-space Ewald and
van der Waals cutoffs 12 Å, electrostatic and polarization Ewald alpha 0.3 Å⁻¹,
32×32×32 PME grid, electrostatic/polarization spline order 5, FFTW, MUTUAL,
and the inherited one-thread / polar-eps 1e-8 settings. Conducting TINFOIL is the
native default: do not add `EWALD-BOUNDARY`, which selects VACUUM in this source.
No response, external-field, dispersion-Ewald, or response-matrix-export record
is present. The ON key differs only by `PME-NORM-EXPORT 1`.

## Exact inventory

Two standalone native programs, OFF then ON, one worker, each calls the original
`initial/getxyz/mechanic/gradient` path exactly once. Thus two analytic energy/
gradient evaluations, zero displaced energies, zero optimizations or extra SCF
calls. Both print full double-precision energy and all 18 gradient components.
Require their complete `harness_ef.dat` bytes to match exactly, including every
energy/force token. OFF must produce no export or application artifacts.

ON additionally performs exactly **54 internal reciprocal applications**: 18
basis columns in the first channel, 18 in the second, then 18 packed applications
with unequal cyclic basis columns. These call actual `grid_uind`, `fftfront`,
stored-qfac multiplication, `fftback`, and `fphi_uind`; local vectors pass through
the same fractional/Cartesian maps as `umutual1`. Input/output normalization is
sqrt(alpha), so the tested operator is normalized R. No model uind/uinp, geometry,
permanent fields, coefficients or external field are modified. Only native PME
scratch qgrid is overwritten. The two model polarization arrays and spline/qfac
arrays are checked unchanged after all applications.

The ON program writes the exact current initialized transform, site identities,
positions, alpha, igrid, first two spline rows and all qfac values immediately
before applications. The collector requires exact equality to the exporter
snapshot from the same gradient call. A mismatch stops admission; matching-grid
assumptions are not substituted for this comparison.

## Numerical checks declared before execution

Reuse frozen SNAPSHOT.py for strict export parsing and ordinary-double norm
replay. Each of five native norms and five qfac summaries must agree within
`1e-12 * max(1, abs(replayed))`. Independently build the small 32768×18 sparse-in-
content scatter matrix from exported weights, apply NumPy's complex FFT with
the native unnormalized backward factor M=32768, and gather with its transpose.
All 54×36 native output components must agree within
`1e-11 * max(1, max_abs_predicted_output)`; record actual maximum error and scale.
This preserves both packed channels and any stored-qfac asymmetry, rather than
assuming equal induced dipoles. Report qfac nonnegativity, conjugate-evenness,
zero mode, channel leakage and nominal U_R. No silent symmetrization or clamping.

These are numerical implementation tolerances, not physical tolerances or an
interval proof of the floating program. A passing result admits only this
initialized export/operator correspondence. Full D stability, trajectory use,
application roundoff budgets and full-box memory remain separate.

## Isolation and costs

Link only the standalone harness against the separately built pinned libtinker
and module files using the same compiler/FFTW/OpenMP recipe. Do not edit its
source tree or replace analyze/testgrad. One immutable output directory per
program; `status='new'` for harness/export artifacts. Stop after any failed
program or missing/incomplete artifact; preserve the other slot as untouched.
Record program attempts, gradient calls and internal reciprocal applications
separately, along with actual wall/child CPU/peak RSS and truthful login identity.
The failed first compile remains authoritative history and is not a native call.
