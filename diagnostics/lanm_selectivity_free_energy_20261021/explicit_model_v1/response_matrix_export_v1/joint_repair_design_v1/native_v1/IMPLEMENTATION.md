# Isolated regularized local-response CPU source

This source implements the reviewed eight-response-amplitude functional in
`finite_screen_v1/interpretation_v1/SPEC.md`, together with the already-existing
native radius/epsilon overrides. It starts from the admitted original-baseline
export source_v2, not the D/M prototype. Original engine files and binaries remain
unchanged. No compiler or molecular program was run during source preparation.

## Input and scope

```
POLAR-LOCAL-RESPONSE 500 36 a_La b_La c_La d_La
POLAR-LOCAL-RESPONSE 501 36 a_Dy b_Dy c_Dy d_Dy
```

The last four entries are physical amplitudes. The evaluator maps root's box
coordinates to c/d using the frozen ln4 transform. The parser rejects duplicate,
missing, nonfinite, out-of-domain and unsupported-pair records. Intrinsic alpha,
pdamp and POLPAIR are unchanged. Absent or all-zero records return before new
mode guards, allocation or arithmetic.

Enabled execution is restricted to one metal or metal-free intact waters, 1–28
positive-alpha active sites, metal monopole exactly +3, vacuum original Thole
MUTUAL without lists, finite cutoffs, prediction or alternative polarization
routes. Static EXTERNAL-FIELD is explicitly supported; frequency-dependent
fields and mixed D/M keys are rejected. Pure water and bare ions have exactly
zero added h/H. This is not a protein, periodic or GPU implementation.

## Five changed files

| File | Change |
| --- | --- |
| polar.f | Opt-in flag, eight coefficients and private geometry arrays |
| induce.f | Parse/guard and rebuild h/H/Jacobians every solve; add h in dfield0a and subtract H*x in ufield0a; native spectral check/readback |
| epolar3.f | Add only −K/2*u·h to analysis polarization energy, per-atom and intermolecular accounting |
| epolar1.f | Add K/2[u·H'v−(u+v)·h'] to analytic gradients and its symmetric virial |
| moments.f | Enabled total moment uses half(uind+uinp); original branch unchanged |

The original reduced energy in epolar0e/epolar1e already contracts corrected
udirp/alpha. It needs no additional H energy. Direct fields retain original
alpha multiplication; all native PCG calls see T*x−H*x through the shared
ufield0a hook. Original positive preconditioning, convergence criterion and peek
remain. No damping scale is inferred by changing polarity. The full geometry
chain includes derivatives through every other water and the metal reaction.

Before each enabled solve, the shared native operator is probed and normalized
by sqrt(original alpha). Existing native `jacobi` supplies eigenvalues/vectors;
symmetry, eigen residual and orthogonality checks reject arithmetic failures
(1e−10 absolute symmetry/orthogonality and scaled eigen residual). The original
SPD criterion remains min>1e−6*max>0. No matrix repair or alternative solver is
used. This private pre-solve probe restores shunt/ncell and zeros trial u/v;
induce0a immediately performs the same normal zero initialization. Prediction
is unsupported. All scalar/domain checks are prospective numerical guards,
not physical-accuracy tolerances.

Only requested fixtures write matrix files. The admitted exporter body is
unchanged: exactly two `response_matrix_001.dat`/`002.dat` snapshots for analyze
EM, including actual At and d/p through the modified shared kernels. The native
export restoration check remains. Testgrad input must omit the export keyword;
its potentially many energy evaluations would exceed that two-snapshot scope.

## Enabled native readback contract

Read one block per induction solve; analyze EM performs two solves. Numeric
values use `es26.17e3`. Whitespace is insignificant. Every expected atom/site is
retained. Pure-water rows have no ENVIRONMENT or PAIR lines; bare ions have a
zero ENVIRONMENT line and no PAIR lines. All-zero/absent keys have no new lines.

```
LOCAL_RESPONSE_PARAMETERS metal_type water_type a b c d
LOCAL_RESPONSE_PAIR metal_atom O_atom gamma0 kappa r_A
LOCAL_RESPONSE_ENVIRONMENT metal_atom s eta sigma beta_Aminus3
LOCAL_RESPONSE_SPECTRUM lambda_min lambda_max symmetry_max eigen_residual_max orthogonality_max
LOCAL_RESPONSE_INDUCED atom u_x u_y u_z v_x v_y v_z
LOCAL_RESPONSE_DUAL max_abs_Ed_minus_Ep max_abs_u_minus_v
```

Fields are in e/Å²; u/v in eÅ. Parameters report all supplied allowed type records;
PAIR reports each actual metal/O pair using native resolved gamma/pdamp. Full
INDUCED vectors permit the conjugate half-sum even for non-export fixtures.
Keep .uind as the original u-only Debye output and an independent formatting
check; do not silently equate it with the dual average. Reconstruct permanent
moments from existing pinned local-frame helpers and preserve common origins.
The native COM-centered M display remains an independent low-precision check.

## Reproducibility and remaining work

The source manifest pins the inherited exporter tree, all candidate source and
CMake files, the five-file patch, specification/domain and helper. Replaying
`patch -p1` from the candidate root onto a copy of the pinned baseline tree
reconstructs the candidate. No new build framework is required: use the already
pinned `metal_pair_response_v1/BUILD.py`/`BUILD.sbatch` with this manifest and a
fresh output directory. They retain `/usr/bin/f95`, RELWITHDEBINFO `-O2 -g
-fbacktrace`, original FFTW and bounded Slurm compiler resources.

The 26 local synthetic/source controls check shape derivatives, water permutation,
rotation, finite contact evaluation, the global coefficient box, exact changed
file set and unchanged exporter text. They are not native numerical qualification.
Planck owns the declared 57-program/636-internal-FD qualification. Root owns
compilation, binding and execution. A successful compile alone does not authorize
fitting. Actual energy/force/field consistency and unaffected-water checks remain
required before the separately declared joint development fit.
