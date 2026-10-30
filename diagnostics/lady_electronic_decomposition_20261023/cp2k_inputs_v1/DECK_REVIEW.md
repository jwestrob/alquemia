# Final quartet deck and MPI source review

Status: **no source-level blocker identified in the four frozen decks**. This is an independent source/metadata review, not molecular execution or numerical admission. No deck, generator, baseline asset or native source was changed. The four-deck generator was read, not rerun. A short standard-library read-only check verified file hashes and literal coordinate/molecule-label tokens.

## Reviewed artifacts

| Artifact | SHA256 |
|---|---|
| `PREPARE.py` | `77efcd979292b4594092969e4e4137afcf1dc7ce87485d4f7a631cda84a9923b` |
| `DECK_PROVENANCE.json` | `916c45531a5384a60dd450d14831e27a6d2efc693eb8d43202a2942e6ada6687` |
| `ENERGY_LEDGER.json` | `24c699cc97c3c55857c56a5ee07dda6213eeaf12e59af69395f23bbd05b732c5` |
| `decks/aquo8_central_La.inp` | `15b07a8ba0a5cb8037847908dd87c6977fa50ed53c8a07ca7489fc7a5d92883b` |
| `decks/aquo8_central_Dy.inp` | `2fde48a9902465724f54dcc8c0b2a3b56fe27cce3ddcd72fd5196e01c3bca46f` |
| `decks/hydrated_acetate_mono_central_La.inp` | `9201271ccceca6d6e8979b9f6e7ea858dfbf3b6ecc75100ec152c2b3c17dd105` |
| `decks/hydrated_acetate_mono_central_Dy.inp` | `33183f873ab07adad9808a4fbe478d265e27a535d1cd23f86bc3f826d0064cb4` |

All three asset hashes match the deck provenance: BASIS_LADY `2ab7190ac0bf52f389fed017672c49f747d75d6d0e92f61dbf4e2f5981b0eb66`; POTENTIAL_LADY `ba87181563f5ca2b8e6e4ff30f43d274752b9c663325905ee6c1dc41c7beccd5`; asset PROVENANCE `ff923ac1a23676e0cd799a989e17271bf30e10aefb74f7a503a50899b82553b4`. The runner must expose those exact files at the literal `/assets/` paths. Metal coefficients are archived literals; H/C/O equivalence remains the documented named-library/BSE correspondence, not a recovery of unprinted historical ORCA coefficients.

## Geometry, state and settings

Every original element and coordinate token is preserved, including negative coordinates and metal at zero. Coordinates are unscaled Angstrom by native default; CENTER_COORDINATES is explicitly false. USER connectivity groups the contiguous **METAL** and **SHELL** names. The trailing numeric `1` and `2` tokens in COORD are ignored by this native format; they are not the authority for grouping. This still yields the intended two fragments, because names differ and each occupies one contiguous atom block (`input_cp2k_subsys.F:977–984`, `topology_generate_util.F:237–275,358–389`). No covalent partition of the shell is introduced.

| Motif, either metal | Atoms | Total charge/electrons | Fragment electrons | Spherical AOs |
|---|---:|---|---|---|
| Aquo8 | 25 = 1+24 | +3 / 88 | metal 8, shell 80 | 391 = 47+344 |
| Acetate + 7 waters | 29 = 1+28 | +2 / 110 | metal 8, shell 102 | 490 = 47+443 |

Every effective calculation is restricted singlet. Fragment charges are +3 and 0/−1. The archived Dy physical-multiplicity metadata is not substituted for its frozen-f effective singlet. The single-token POTENTIAL names correctly select ECP from their prefix (`qs_kind_types.F:1868–1888`); semilocal ECP is explicitly enabled, with La46/Dy55 and the already reviewed guarded Dy occupation reader.

The actual decks match the fixed settings reviewed in PILOT_SETTINGS_REVIEW.md: 30-A nonperiodic MT cell; GAPW 1200/100-Ry five-level grid; 200 radial/590 angular points; explicit tight Gaussian/GAPW thresholds; ordinary fragment OT at EPS_SCF 1e-10; molecular fresh ALMO guess, DIAG block optimization and FULL_SCF/PCG release at EPS_ERROR 1e-9; no extrapolation, truncated virtual space, ADMM or finite-range exchange. PBE exchange/correlation scaling 0.75/1 plus full Coulomb HF fraction 0.25 is PBE0. D4 is absent from the CP2K electronic calculation and belongs only in the explicit unchanged archived-energy ledger.

The root's fixed gates govern: **1e-4 Eh cross-code full-energy agreement and 1e-8 Eh arithmetic/order checks**. The earlier illustrative 0.1-kcal value is not the execution gate. Neither the settings nor output closure guarantee these gates or physical accuracy.

## MPI distribution and output checks

Actual release source is the already pinned archive under `workspaces/lady_electronic_decomposition_20261023/cp2k_build_v1/cp2k-2026.2`, revision git:c92cc08. `mscfg_methods.F:186,239–306` obtains the parent parallel environment and runs a serial loop over fragments. Each `create_small_subsys` and `qs_init` receives the **same full MPI communicator**. Thus 12 MPI × 2 OpenMP uses the 12 ranks for each fragment in turn; it does not allocate one rank per fragment or split ranks 6/6. There is no requirement that 47, 344, 443 or the total AO counts divide evenly by 12. `mscfg_types.F:123–180` assembles variable-size fragment matrices using their actual row/column sizes and checks the total dimensions. Tiny-metal work can scale poorly; no throughput benefit is established by this source review.

`mscfg_methods.F:397–410` prints the unambiguous `MOLECULAR GUESS: FRAGMENT i OUT OF n` banner, charge/multiplicity, and original atom-index/name list. Initialization metadata can occur before that banner, because qs_init precedes it. Each following fragment SCF emits electron/occupied-orbital/AO counts, the explicit converged marker, and its `Total energy:` through PROGRAM_RUN_INFO, whose default print level is LOW and is active at MEDIUM (`input_cp2k_scf.F:351–355`; `qs_scf_initialization.F:1179–1183`; `qs_scf_output.F:143–178,873`; `qs_scf_loop_utils.F:594–601`). The two fragment segments must be kept separate from the parent ALMO DIIS/PCG records. Reject missing/failed convergence independently of telescoping energy arithmetic.

Retain fragment energies, their printed sum R, initial antisymmetrized energy F, final converged DIIS energy P, and final full PCG/ENERGY U. The identities are FRZ=F−R, POL=P−F, CT=U−P. FULL_SCF does not require the perturbative helper's `ENERGY OF BLOCK-DIAGONAL ALMOs` label. Only combined frozen interaction is supported by this route; no separate electrostatic/Pauli split is inferred. Geometry/state/counts, convergence, closure, monotone variational ordering and full-energy agreement remain actual-run admissions.
