# Donor source handoff

The isolated tree starts from the exact admitted water source manifest82886b76. `PREPARE_SOURCE.py` verifies its file pins, copies source/CMake into a fresh versioned workspace, and changes only **polar.f, induce.f and testgrad.f**. `DONOR_RESPONSE.patch` reconstructs those changes. No original source, executable, donor key or parameter file was edited. No build or molecular program has run for this version.

## Integration

- `polar.f` adds a donor-active flag, four separate coefficients and the oxygen beta derivative tensor. The water fields/coefficients and original alpha/pdamp remain unchanged.
- `localresponsesetup` resets the original state, calls the donor setup, and returns early only for a matched nonzero donor channel. All original water parsing, geometry, spectrum and arithmetic remain intact. Absent/zero/unmatched donor records return before donor mode guards. Nonzero water response therefore remains available when no O403 exists, even if donor records are in the same key.
- The donor graph validator requires the exact acetate401/402/403/404 connectivity, treating both oxygens identically. The enabled channel requires one unbonded500/501 metal and eight real active polarizable sites. The exporter additionally recognizes the same exact seven-site ligand. No mixed ligand topology is admitted.
- Existing `dfield0a` adds h to both fields; existing `ufield0a` subtracts beta*x at each corrected site. These routines are byte-identical to the parent tree. The native spectrum, PCG and two-snapshot exporter body are unchanged.
- The new donor gradient sums the beta contractions over both oxygen sites, adds h derivatives, and includes their opposite metal derivatives and symmetric virial. Dispatch leaves the old water gradient arithmetic unchanged. `epolar.f`, `epolar1.f`, `epolar3.f`, `moments.f`, `rotpole.f` and `torque.f` are byte-identical to the parent. Their corrected-u/v energy and moving-frame derivative paths are reused.

## Keyword and readback

```
POLAR-DONOR-RESPONSE metal_type 403 a c
```

Native accepts finite amplitudes with arithmetic/topology/state/SPD guards. The qualification runner, not the engine, restricts executed coefficients to the six prospectively declared drive/stiffness/combined0.1 states. Extra, malformed or duplicate records fail. A donor channel actually active on acetate rejects the other response-family keywords; a no-donor bypass does not reject the retained water channel.

One set is emitted per induction solve; analyze EM produces two solves. Real values use `es26.17e3`:

```
DONOR_RESPONSE_PARAMETERS metal_type 403 a c
DONOR_RESPONSE_PAIR metal_atom O_atom gamma kappa r_A f w beta_Aminus3
DONOR_RESPONSE_FIELD O_atom hx hy hz
DONOR_RESPONSE_PERMANENT_LAB_EA mux muy muz
DONOR_RESPONSE_CHARGE total_charge
LOCAL_RESPONSE_SPECTRUM min max symmetry_max eigen_residual_max orthogonality_max
LOCAL_RESPONSE_INDUCED atom ux uy uz vx vy vz
LOCAL_RESPONSE_DUAL max_abs_Ed_minus_Ep max_abs_u_minus_v
```

h is in e/Å²; u/v and the permanent laboratory dipole are in eÅ. Permanent reporting occurs inside induce with freshly rotated permanent rpole, **before** moments adds induced dipoles. It uses sum(q*r+rotated permanent dipole) at fixed laboratory origin0, with no COM subtraction. It is also emitted for exact donor/ion exports when the correction is absent or zero. Existing exporter file format and `.uind` u-only output remain unchanged. Compute the conjugate dipole from permanent+half(sum(u+v)); do not assume donor u=v or interpret the zero native M display as conjugate.

Immediately after testgrad's complete analytic gradient, an active donor emits:

```
DONOR_RESPONSE_TOTAL_ENERGY E_kcal_mol
DONOR_RESPONSE_TOTAL_GRADIENT atom gx gy gz
DONOR_RESPONSE_TOTAL_VIRIAL row v1 v2 v3
```

These are the full native gradient (kcal/mol/Å) and virial (kcal/mol), including original bonded/nonbonded and moving-frame terms. They are read-only observations, not correction-only replacements. Zero/off testgrad retains the original output path. Testgrad must omit RESPONSE-MATRIX-EXPORT because its energy differences exceed the exporter two-snapshot scope.

## Checks, build and remaining evidence

SOURCE_CONTROLS verifies pinned sources, unchanged water/operator/exporter/frame paths and synthetic unequal-dual-field derivatives, oxygen permutation, translation/reactions, contact finiteness and drive/stiffness susceptibility structure. Synthetic arrays are algebra checks, not molecular qualification or target evaluation. No compiler was invoked.

READY pins the existing `metal_pair_response_v1/BUILD.py` and `BUILD.sbatch`, with this new manifest/hash and a fresh output directory. Root owns the twelve-CPU Slurm compile. Expected products are `build/analyze`, `build/testgrad`, `build/libtinker.a` and build receipts. The manifest’s `source_outputs`/`cmake` interface is unchanged, so no new build framework is needed.

Independent source review precedes compilation. Arendt owns the separately declared121-program/576-internal-Cartesian-difference numerical matrix. No molecular execution, fit, new physical threshold or coefficient selection is released by this source handoff. Failed Dy field qualification and pending reference-method evidence remain explicit prerequisites to scientific interpretation and fitting.
