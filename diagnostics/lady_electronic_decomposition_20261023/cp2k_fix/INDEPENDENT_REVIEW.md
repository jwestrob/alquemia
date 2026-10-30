# Independent source review: Dy55 occupation patch

2026-10-23. **Passed source review for the bounded archived Dy55 use; no correction requested.** No installation, build, electronic calculation or scientific-runtime import was performed. The only executed check was a standard-library, in-memory text application/hash comparison; the parent source was not edited.

## Exact reviewed inputs

CP2K v2026.2, commit `67b5da876dd6a76b8b021d5a04d1c81ba79a4c50`:

| Artifact | SHA256 |
|---|---|
| Parent `src/atom_types.F` | `1a91f3c5da4403c62db24d37cad13aa47873459d61d9f64899697e82121c3ee9` |
| `dy55_default_occupation.patch` | `7a4fcb3d977f90b108bb1f1f61b45c78e4a2e4236d18b122df9b9f9e96f265e2` |
| Prospective patched `src/atom_types.F` | `9adc98b30524851370c17aad8f4a9b6a95e1ca8b0ae50c21c479029d2e55c60a` |
| Inspected `src/qs_kind_types.F` | `40ec254a7820bb3976ff7661db1053d11af8a5d33a231f9123a81c5e210736b1` |
| Inspected `src/common/periodic_table.F` | `b8971d522856d6b6e4eaba2721684cfbdde46474fb5b83f0556cbf7f90ecefa8` |

The patch has exactly seven added lines and zero removed lines. Its context matches the parent once; the resulting hash agrees with the author's prospective hash.

## Scientific and implementation finding

In this version's `read_ecp_potential`, `ncore` is the atomic number obtained from `get_ptable_info`; `nel` is the removed-electron count read from the ECP. Thus `CASE (55)` with `ncore /= 66` rejection is correctly restricted to Dy55. Existing effective charge remains `zion = 66 - 55 = 11`.

The new default subtracts `[Kr] 4d10 4f9` from neutral Dy's angular occupation totals `(12,24,20,10)`, yielding `(4,6,0,1)`, eleven explicit electrons. Molecular charge +3 independently gives the required eight-electron metal fragment. The added branch does not redefine the molecular singlet or insert a physical valence-f claim. [Exact parent reader](https://github.com/cp2k/cp2k/blob/67b5da876dd6a76b8b021d5a04d1c81ba79a4c50/src/atom_types.F), [Dy table entry](https://github.com/cp2k/cp2k/blob/67b5da876dd6a76b8b021d5a04d1c81ba79a4c50/src/common/periodic_table.F)

All local/projector coefficients, radial powers, effective charge, existing core cases and the final nonnegative-occupation assertion remain unchanged. The patch applies after either supported ECP input-reading path. Downstream potential construction uses the original arrays, and `set_pseudo_state` already supports lanthanide cores 47–60 with `4f = Ncore - 46`; Dy55 therefore reaches the existing nine-core-f-electron path. [Exact downstream source](https://github.com/cp2k/cp2k/blob/67b5da876dd6a76b8b021d5a04d1c81ba79a4c50/src/qs_kind_types.F)

This preserves the intended formal semilocal-ECP Hamiltonian. Occupation bookkeeping still affects atomic initialization and can affect convergence or the attained SCF solution. Source equivalence therefore does not establish full-energy agreement, ALMO convergence, or decomposition closure. Root must bind the actual patched build and unchanged translated inputs before molecular execution; the existing quartet gates remain necessary. This review establishes neither general core55 support nor a qualified La/Dy result.
