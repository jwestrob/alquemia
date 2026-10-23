# Original-baseline response-matrix export

The saved oxygen-only Jacobian is locally too weak for the four observed state
curvature differences in the declared domain. Before implementing a correction,
export the **original** native response operator and both driving fields. This
permits exact finite linear-model calculations and stability assessment on
consumed geometries without another quantum panel or a speculative force patch.
It does not validate a proposed correction or supply coordinate derivatives.

## Minimal isolated source change

`PREPARE_SOURCE.py` copies the original `tinker_source/source` and `cmake` into
an immutable `response_matrix_export_v1/source_v1`. Only `induce.f` changes:

1. At `induce` entry, recognize the no-argument `RESPONSE-MATRIX-EXPORT` reporting
   keyword and check the supported native route. Absent keyword returns before
   any guard or numerical work. No global electrostatic flag is changed.
2. Immediately after the original induction solver returns, before predictor
   history and native moment printing, export the operator and restore state.
3. Append two reporting helpers in the same compilation unit. Original
   `induce0a`, `dfield0a`, `ufield0a`, PCG, damping, energies and gradients remain
   byte-for-byte unchanged. The D/M prototype is not the source or binary.

The exporter accepts original nonperiodic, no-list, original-Thole MUTUAL PCG,
no solvent/CHGPEN/EXCHPOL/charge flux/predictor/oscillating field, and no finite
multipole cutoff. All current sites must be active, have alpha≥1e−8 Å³ (native
PCG `polmin`), and be types36/37/500/501, at most one metal. n=npole is between
1 and28; bare ions are admitted. Static external fields remain usable, with
their actual native `texfld` recorded. With no reporting keyword the original
field path and all other original modes remain unchanged.

## Exact operator, right-hand sides and restoration

Native `ufield0a(field,fieldp)` returns T*uind and T*uinp. The original PCG loop
forms `conj/poli−field`; its unpreconditioned physical operator is therefore
`A=diag(1/alpha)−T`, not its preconditioner. The alpha guard ensures the native
`max(polmin,alpha)` equals alpha in every exported row.

Save all uind, uinp, udir and udirp values. For each Cartesian unit vector in
the exact `ipole` order, set both temporary induced arrays to that vector and
call the **existing** `ufield0a`. Negate its first channel, add the diagonal
inverse alpha, and record the channel difference. No damping formula or
exclusion table is reimplemented. A is not silently symmetrized.

Export the actual solver fields `Ed=udir/alpha`, `Ep=udirp/alpha`, and converged
u/v. Independently recompute the direct fields through original `dfield0a` plus
the recorded applied field, for a closure diagnostic. After restoring u/v,
apply the native mutual operator once more and compare its action to dense A*u
and A*v. Retain Au−Ed and Av−Ep separately; a finite solver residual is not
assumed to bound solution error without matrix conditioning.

`dfield0a/ufield0a` write local arrays and invoke `switch('MPOLE')`. This changes
the 18 shunt values and calls `replica`, which sets `ncell=1` before returning
for a nonperiodic system. The exporter snapshots/restores those values and
ncell as well as all four induced/direct arrays, and checks exact restoration.
The guarded kernels do not modify permanent multipoles, geometry, exclusion
tables or predictor histories. Predictor/periodic/replica routes are rejected
before any solver call when export is requested.

One analyze-EM program makes two native induction solves. Export both as
`response_matrix_001.dat` and `response_matrix_002.dat`, each opened with
`status='new'`; never overwrite. A third solve is rejected in this diagnostic.
Basis operator applications are internal reporting operations, not molecular
energy calls or Cartesian finite differences.

## Public text contract

All floating fields use `ES26.17E3`; ordering is explicit:

```
TINKER_RESPONSE_MATRIX_V1
COUNTS sequence n npole ndim
UNITS A=Angstrom^-3 RHS=e/Angstrom^2 U=eA
CONSTANTS electric dielec debye
EXTERNAL use_exfld texfld_x texfld_y texfld_z
CHECKS dual_operator_max direct_RHS_max symmetry_max operator_action_max
       Au_minus_Ed_max Av_minus_Ep_max       # all six on one record
RESTORED T
SITE ipole_position atom_index type native_name alpha_A3 x_A y_A z_A
... npole SITE records
VECTOR dof Ed Ep u v
A row A_row_1 ... A_row_ndim
... ndim VECTOR/A pairs
END_TINKER_RESPONSE_MATRIX_V1
```

The exact site map, coordinates and alpha accompany every file. A and both
fields use native charge/Å units; energy is `−(electric/dielec)/2 * uᵀEp`.
For a static field, Ed/Ep include texfld; subtract its repeated site vector
only when reconstructing the field-free driving field. Uniform susceptibility
uses the energy-consistent field conversion already pinned by model19, not a
new unit constant. Dense inverse/matrix solving and stability analysis belong
to the saved collector, not to this reporting hook.

## Build and numerical qualification boundary

Reuse the pinned existing `metal_pair_response_v1/BUILD.py` and `BUILD.sbatch`
with this source manifest. Their existing CMake recipe uses the copied `cmake`,
RELWITHDEBINFO, `/usr/bin/f95`, `-O2 -g -fbacktrace`, system static FFTW and
OpenMM OFF; build analyze/testgrad in a fresh root-owned directory. No new
build framework or original binary replacement is needed.

First scout: the exact consumed W9 La/water01 dimer used by model19, with the
same selected17676 key/parameters and coordinates. Two new analyze-EM programs
are proposed: rebuilt binary with export absent, then with the single reporting
keyword. Compare both to admitted zero-field energy/components and saved site
dipoles; use the four admitted ±h/±h/2 fields without rerunning them. Planck owns
the small scout adapter/declaration and its prospective numerical tolerances.
It must verify complete two-file denominators, restoration, matrix symmetry and
channel equality, native/dense action, solve and polarization-energy replay,
condition-aware residual bounds, and the existing uniform-z response column.
No forced symmetrization, unsupported-mode approximation or tolerance tuning.

The exact larger inventory is Singer's separate preparation (currently54
consumed states, including2 bare ions); root declares/releases it only after
the scout audit. No finite h/H application, fit or native correction is part of
the exporter. Exported A/Ed/Ep do not provide dA/dR or dEd/dR, and cannot
establish analytic forces for a future functional. No production, GPU,
periodic, donor or held-out release follows.

## Preparation status

Source-copy assertions verify that only induce.f changes, helper code fits
72-column fixed form, and every original solver/kernel remains intact.
MANIFEST.json pins all original and output source files. This is uncompiled
source preparation, zero native calls. Root owns compilation and execution;
actual numerical parity is pending the independent scout.
