# Numerical La-to-Dy parameter path

Declared 21 October 2026 before preparation or execution. This is a numerical
implementation test using the preserved, physically unqualified candidate
parameters. It does not release production sampling or fit any parameter.

## Question

Can separate native parameter files represent a continuous La-to-Dy path with
exact endpoint energies and forces, including pair-specific van der Waals and
polarization terms? This avoids assuming that the engine's decoupling keywords
already implement substitution between two interacting metal species.

Use the two existing, CPU/GPU-checked fixed fixtures: water_central and
acetate_bi_central. Retain every coordinate, bond, nonmetal type and nonmetal
parameter. The existing La500 and Dy501 states are the endpoints. Introduce one
artificial type/class502 with fixed charge +3, zero higher permanent multipoles,
fixed mass138.91, and no bonded neighbors. Its atomic number57 is a bookkeeping
label, not a fractional-element quantum calculation.

The metal mass is held fixed throughout this configurational path. The omitted
metal kinetic mass term cancels between bound and aqueous exchange legs when
the same classical unconstrained degrees of freedom and mass convention are
used in both. This does not cover isotope/nuclear-quantum effects or mismatched
mass-dependent constraints. This test contains no dynamics or free energy estimate.

## Path definition

For lambda0=La and lambda1=Dy, linearly interpolate metal polarizability and
each resolved metal/nonmetal pair's radius, well depth and Thole damping.
Resolve endpoint pair overrides first; where absent, use the exact pinned
native CUBIC-MEAN radius, HHG well-depth and minimum-nonzero Thole rules.
Explicitly write each metal/nonmetal class VDWPR and type POLPAIR in the new
key. Do not interpolate an existing La override toward zero, drop it abruptly,
or assume that a type identifier and a class identifier are interchangeable.

Native source authority is the successful GPU build's pinned Tinker submodule
`kvdw.f` (diameter-to-radius conversion, CUBIC-MEAN, HHG, VDWPR semantics) and
`kpolar.f` (default and pair-specific Thole damping). Require the expected
radius/epsilon conventions; reject unsupported recipes rather than guess.

The path is a defined interpolation of force-field parameters, not a linear
interpolation of minimized endpoint energies. Mutual induced dipoles must be
solved anew at every evaluation. Positive radii/well depths/polarizabilities and
constant +3 charge avoid introducing a disappearing-particle endpoint here.

## Finite numerical inventory and acceptance

For each fixture prepare lambda0,.25,.5,.75,1 and .5 plus/minus1e-4 and5e-5:
18 fixed GPU energy/gradient evaluations total. At the two midpoints also
evaluate the unchanged CPU engine, for20 native program invocations total.
Use the existing double-precision H200 build/shared CUDA12 runtime. Root owns
submissions. No native calls occur during source/preparation validation.

- Endpoint energies agree with the saved corresponding native endpoints to
  1e-4 kcal/mol; maximum force-component discrepancy <=1e-3 kcal/mol/angstrom.
- Midpoint CPU/GPU energies and forces meet the same tolerances.
- Spatial finite-difference maximum error <=1e-3 kcal/mol/angstrom at the two
  midpoints, using the already qualified displacement convention.
- Central lambda energy derivatives at the two step sizes agree within
  1e-3 kcal/mol. Both are finite-difference estimates; do not label either an
  analytic derivative or use this small fixture test as a convergence proof.
- Every stage has finite arrays, exact atom/state accounting, actual native
  completion and source pins. A native zero exit without completed output is
  failure. Retain all failed attempts and do not loosen gates after inspection.

Outputs are a reusable parameter renderer, finite test record and explicit
limitations. Periodic protein/water electrostatics, dynamics stability, ensemble
sampling, overlap and physical-model qualification remain separate work.
