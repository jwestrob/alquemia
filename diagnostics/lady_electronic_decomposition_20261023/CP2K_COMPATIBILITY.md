# CP2K compatibility with the exact ALMO quartet

Source/documentation inspection only, 23 October 2026 UTC. No installation, source patch, electronic calculation or scheduler action was performed.

## Finding

**CP2K is a genuine free native-ALMO implementation, but the inspected stock versions cannot read the quartet's original Dy 55-core ECP through their normal ECP input route.** This is a concrete parser/default-occupation restriction. It is not evidence that semilocal ECP integrals require a different Hamiltonian.

The same parser accepts La46. The PBE0/ALMO/semilocal-ECP ingredients have a plausible common execution path, but no combined La or Dy calculation has been tested here. Installing stock CP2K alone would not remove the identified Dy55 blocker.

Local executable and entitlement inspection is recorded separately in [SOFTWARE_ACCESS.txt](SOFTWARE_ACCESS.txt) and [SOFTWARE_ACCESS.json](SOFTWARE_ACCESS.json). No CP2K executable was found within that bounded inspection; this is not a campus-wide absence claim.

## Exact source versions

- Stable `v2026.2`: commit `67b5da876dd6a76b8b021d5a04d1c81ba79a4c50`.
- Inspected master: commit `62906da87cd0e92000c341d4500e1ed4a5a6af9c`, dated 2026-10-23 18:00 UTC.
- `v2026.1` and `v2025.2` were also checked specifically for the core-count restriction.

Public-source hashes are retained in `SOFTWARE_ACCESS.json`. These reads do not select a build version or authorize a pilot.

## Hard blocker: original Dy55 input

In stable [atom_types.F, lines 3033–3068](https://github.com/cp2k/cp2k/blob/67b5da876dd6a76b8b021d5a04d1c81ba79a4c50/src/atom_types.F#L3033), `read_ecp_potential` constructs default atomic occupations using a core-electron-count switch. The allowed counts are:

`0, 2, 10, 18, 28, 36, 46, 54, 60, 68, 78`.

Core55 reaches `Unknown Core State`. Stable `v2026.2` and `v2026.1` have identical `atom_types.F` SHA256 `1a91f3c5da4403c62db24d37cad13aa47873459d61d9f64899697e82121c3ee9`; `v2025.2` has the same restriction. Master retains it at lines 3047–3082.

Both inline and file-based ECP input converge before this unconditional switch. `KIND/ELEC_CONF` does not bypass it: the ECP reader must first return, after which the kind's explicit/default occupations are applied. See master [qs_kind_types.F, lines 2799–2871](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/qs_kind_types.F#L2799). No documented exact-input alternate route was found.

Changing the declared core to54 or60, substituting an element or GTH potential, or applying an unverified charge correction would not be an admitted replay of the original input. None was proposed or performed. A separately authorized implementation change would require its own review; this assessment contains no patch.

## What the source does support

### Original semilocal potentials and Gaussian bases

The official [KIND documentation](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/SUBSYS/KIND.html) supports Gaussian orbital bases and `POTENTIAL ECP` with `ECP_SEMI_LOCAL TRUE`. This preserves the original semilocal form rather than forcing a GTH replacement.

The archived metal assets contain scalar Gaussian semilocal potentials, s–f nonlocal channels, a zero g local channel, and contracted s–f orbital bases:

- La asset SHA256 `c5298b309b2d897e5de59ddb6b4c7992913dd8cba18b7cd21f9757ce09b167d3`.
- Dy asset SHA256 `87bb209a42eca365c93d8f69088a57371285623a1d4afa1c1336390279fafb20`.
- Both are under `workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/common_basis_decomposition_v1/prepared_v1/assets/`.

These data need a literal format conversion preserving effective charge, coefficients, exponents, angular channels, radial-power conventions and orbital normalization. Zero local coefficients are not rejected by the inspected read loop. This does not overcome its later core55 rejection. Preserve both effective-singlet M3+ fragments with eight explicit electrons; do not substitute a physical Dy sextet calculation.

### PBE0 and native ALMO

CP2K documents [PBE0 exact exchange and native ALMO](https://www.cp2k.org/features). More concretely, master [almo_scf_qs.F, lines 639–680](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/almo_scf_qs.F#L639) sends the ALMO density through the common Quickstep KS update. That update adds exact exchange via `hfx_ks_matrix` when requested; see [qs_ks_methods.F, lines 1085–1119](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/qs_ks_methods.F#L1085).

No categorical ALMO/HFX prohibition was found in the inspected path. This is source support, not a combined PBE0/custom-ECP qualification. The inspected ALMO regression examples use semilocal functionals/GTH; separate ECP examples do not establish the combined quartet.

The [ALMO_SCF documentation](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/ALMO_SCF.html) distinguishes the constrained block-local state from `FULL_SCF`, which restores unrestricted spatial delocalization. Its default `XALMO_X` is a single-excitation correction with spatial restrictions and is not equivalent to the requested full relaxation. A total CT energy would be the admitted full-state minus constrained-state energy; an approximate pairwise CT table must not silently replace that difference.

### Frozen energy without trusting a charged-periodic summary

The [ANALYSIS documentation](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/ALMO_SCF/ANALYSIS.html) explicitly warns that its `SUBLATTICE` frozen-reference term lacks interpretation for charged fragments. Its wording does not limit this warning to periodic calculations. The exact quartet contains charged fragments, so it cannot be ignored.

The raw-energy source path is more informative:

1. [mscfg_methods.F, lines 174–182 and 291–319](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/mscfg_methods.F#L291) reuses the parent DFT/FORCE_EVAL settings and parent cell for each fragment, with explicit fragment charge/multiplicity. It stores each fragment's total energy. The code does not force periodicity; removing periodicity automatically is commented as unfinished work, so the parent boundary conditions must already be correct.
2. [almo_scf.F, lines 626–638](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/almo_scf.F#L626) assembles those fragment orbitals. Lines747–788 construct the nonorthogonal antisymmetrized projector and evaluate its energy before ALMO relaxation.
3. [Lines790–795](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/almo_scf.F#L790) print the summed fragment energy and initial-guess energy independently of the formal `ANALYSIS/SUBLATTICE` summary.

Thus a nonperiodic same-cell calculation has native raw ingredients for fragment, frozen, constrained-polarized and full-state energies. This source finding prevents a categorical rejection based solely on the SUBLATTICE warning. It does not prove numerical reference consistency, convergence, or finite-cell independence. The declared fragmentation must be used exactly, including the whole water shell as one fragment where specified.

## Build requirements, if root later chooses a compatible route

- **Semilocal integrals:** LIBGRPP is already bundled in stable [src/CMakeLists.txt, lines1461–1501](https://github.com/cp2k/cp2k/blob/67b5da876dd6a76b8b021d5a04d1c81ba79a4c50/src/CMakeLists.txt#L1461). It is not a separately missing proprietary license. Preserve `ECP_SEMI_LOCAL TRUE`.
- **Exact exchange:** enable `CP2K_USE_LIBINT2` with Libint2's Fortran interface. The stable [CMake checks](https://github.com/cp2k/cp2k/blob/67b5da876dd6a76b8b021d5a04d1c81ba79a4c50/CMakeLists.txt#L817) explicitly require `libint_f.mod`. Verify sufficient orbital angular-momentum coverage for the actual basis.
- **Base runtime:** compatible C/C++/Fortran compilers, BLAS/LAPACK, DBCSR, OpenMP and FFT support; MPI/ScaLAPACK for an MPI build. Native ALMO sources are compiled as part of CP2K; no separate ALMO license was identified.
- **Physics:** preserve vacuum boundaries, exact PBE0 exchange fraction, original basis/ECPs and fragment states. Do not silently enable ADMM, truncate delocalization, or include D4 in an electronic CT contribution. Numerical grids/Poisson settings differ from ORCA and need explicit review rather than an exact-replay claim.

**Decision boundary:** CP2K remains a real free implementation worth considering, but a stock installation is presently insufficient for the unchanged La/Dy quartet because of the confirmed Dy55 reader restriction. No claim about all free implementations or institutional access follows from this bounded inspection.
