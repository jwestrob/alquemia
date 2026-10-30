# CP2K La-aquo ALMO pilot: independent source review

2026-10-23. Source/documentation inspection only; no installation, runtime import, calculation or input execution. Inspected CP2K source commit: `62906da87cd0e92000c341d4500e1ed4a5a6af9c`. Root owns the eventual executable/version and numerical settings.

**Verdict:** the native molecular-guess route supports the required isolated-fragment, frozen, block-polarized and fully delocalized energy chain for a genuinely nonperiodic La pilot. It does not require interpreting the charged periodic `SUBLATTICE` analysis. This resolves the specific reference concern in `METHOD_REVIEW.md`; it is not a numerical qualification. **A separate stock-source blocker presently prevents promising the exact quartet: the ECP reader rejects Dy's required 55-electron core.** Do not spend a La pilot as though that quartet compatibility were already established.

## Why the native energy reconstruction is legitimate

`mscfg_methods.F` retains the parent DFT/FORCE_EVAL settings, creates each fragment from its own atoms in the parent cell, passes explicit fragment charge and multiplicity into `qs_init`, runs `qs_energies`, and stores `energy_of_frag`. Thus a nonperiodic parent also gives nonperiodic isolated fragments at their unchanged complex coordinates, without the other fragment's basis. This is the actual native molecular initialization, not an independently constructed polarization method. [Exact fragment source](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/mscfg_methods.F)

`almo_scf_initial_guess` assembles those molecular orbitals, normalizes/orthogonalizes within their local blocks, builds the nonorthogonal antisymmetrized density through the full overlap matrix, and evaluates its energy. It prints both the summed single-molecule energy and initial-guess energy independently of the optional `ANALYSIS/SUBLATTICE` branch. These are respectively the fragment reference and frozen state when the guess is genuinely fresh and molecular. Atomic guesses or extrapolated/restarted guesses would not establish this frozen reference. [Exact ALMO initialization](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/almo_scf.F)

The block optimizer then supplies the no-metal–shell-delocalization state. `FULL_SCF` releases spatial restrictions and optimizes the full state; default `XALMO_X` is an approximate delocalization correction and is not acceptable here. The optimizer retains the block energy and reports energy lowering/corrected energy; the final native total and convergence record must also be retained. [ALMO controls](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/ALMO_SCF.html), [exact optimizer](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/almo_scf_optimizer.F)

The ALMO density uses the ordinary KS update, whose source includes Hartree–Fock exchange. No exclusion of hybrid exchange was found in this inspected call path. This supports PBE0 compatibility but does not prove the combined ALMO/GAPW/semilocal-ECP/HFX implementation numerically. [ALMO KS interface](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/almo_scf_qs.F), [KS/HFX implementation](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/qs_ks_methods.F)

## Required La-aquo native input contract

Use the pinned 25-atom `inputs/xyz/aquo8_central_La.xyz`, SHA256 `371533e94b5e7b2bbe61f1ecec77d1ddbfdbc4edd0beea98af77dad23bf39c4a`. Preserve all coordinate tokens and ordering. Input-bundle authority is `inputs/PROVENANCE.json`, SHA256 `677c759c1045e2d7d5ab2dd05977e4a516ca04ff8d520c40111669eabb452d66`.

| Setting | Required choice |
|---|---|
| Task/state | Single-point energy; total charge +3, effective singlet, 88 explicit electrons; no smearing or excited-state substitution |
| Domain 1 | Atom 1, La3+, charge +3, multiplicity 1, eight explicit electrons |
| Domain 2 | Atoms 2–25 together, entire eight-water shell, charge 0, multiplicity 1, 80 electrons |
| Grouping | `TOPOLOGY/CONNECTIVITY USER`; coordinate suffixes identify metal group 1 and the complete shell group 2. No automatic bond-based water fragmentation |
| Electrostatics | Both `SUBSYS/CELL/PERIODIC NONE` and `DFT/POISSON/PERIODIC NONE`; isolated Poisson solver such as `POISSON_SOLVER MT`, with an explicitly sufficient cell |
| Native method | `QS/METHOD GAPW`, `QS/ALMO_SCF T`; exact PBE0: 0.75 PBE exchange + full PBE correlation + 0.25 full-range HF exchange |
| ALMO | `ALMO_SCF_GUESS MOLECULAR`, fresh start, `ALMO_EXTRAPOLATION_ORDER 0`, `DELOCALIZE_METHOD FULL_SCF`; explicit `FRAGMENT 1` and `FRAGMENT 2` charges/multiplicities |
| ECP/basis | Exact archived La orbital exponents/contractions and 46-electron semilocal ECP; `KIND/ECP_SEMI_LOCAL T`. O/H use exact all-electron def2-TZVP functions; no GTH substitution, ADMM surrogate or changed contraction |
| Analysis reference | Do not rely on charged `FROZEN_MO_ENERGY_TERM SUBLATTICE`; retain native absolute molecular/frozen/block/full energies |

Explicit molecular grouping and charged singlets are demonstrated by the [official ALMO regression input](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/tests/QS/regtest-almo-1/almo-singlet-charge.inp). The [fragment interface](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/ALMO_SCF/FRAGMENT.html) supplies charge and multiplicity separately. All waters sharing a domain allows their internal delocalization, as the requested whole-shell partition requires.

`ECP_SEMI_LOCAL` means evaluation of the original semilocal potential; disabling it changes its representation. The archived ORCA ECP's radial powers and local/projector channels must be translated literally, not selected by a similar library name. [KIND documentation](https://manual.cp2k.org/cp2k-2024_1-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/KIND.html)

**This is an input contract, not a run-ready deck.** Cell size, grids/GAPW accuracy, integral thresholds and all fragment/block/full stopping criteria remain to be bound prospectively. MT requires a cell exceeding twice the cluster extent including density; default or tiny regression-test cells are not evidence of convergence. Preserve original Cartesian coordinates while choosing the enclosing cell. [Poisson documentation](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/POISSON.html) Cross-program full-energy agreement must qualify the numerical realization of archived PBE0/RIJCOSX, not assume identical implementation.

## Closure and bookkeeping

Let `R = sum(native isolated-fragment energies)`, `F = native molecular initial-guess energy`, `P = converged block energy`, and `U = converged FULL_SCF energy`. For electronic energies:

```text
FRZ = F - R; POL = P - F; CT = U - P
FRZ + POL + CT = U - R
```

One native full run can retain this chain; a separate `NONE` calculation is not logically necessary merely to obtain the already retained block energy. All internal fragment/block/full SCFs nevertheless count as actual work. Require their convergence and the nested variational ordering within declared numerical resolution, not an energy clamp. Printed energy closure is only arithmetic consistency; it cannot certify SCF accuracy.

The archived D4 correction must appear once in the full-energy comparison. Native D3 or an unrelated dispersion decomposition cannot replace it. With the same geometry-only archived D4 constant added to F/P/U, POL and CT are unchanged; the fragment-reference D4 ledger must be specified separately for interaction closure. Preserve electronic and D4 quantities separately. Use uncorrected energies for the raw closure; do not combine a BSSE-corrected CT with an uncorrected total. Own-fragment-basis energies do not contain deformation to separately relaxed fragments.

## Exact quartet blocker and limits

At the inspected commit, `read_ecp_potential_file` recognizes core counts 0, 2, 10, 18, 28, 36, 46, 54, 60, 68 and 78. The default branch aborts with an unknown-core-state error. La46 is recognized; archived Dy55 is not. A later electronic-configuration override does not bypass this earlier reader failure. **Do not replace Dy55 by 54/60 or claim a completed compatible quartet.** Singer owns any documented-path resolution; no source patch or alternative Hamiltonian is proposed here. [Exact ECP reader](https://github.com/cp2k/cp2k/blob/62906da87cd0e92000c341d4500e1ed4a5a6af9c/src/atom_types.F)

If eventually qualified, this would be finite-basis AO-block ALMO: CT is the energetic effect of releasing the metal–whole-shell orbital constraint. It is not a unique microscopic partition, explicit-f analysis, or protein/free-energy inference. The original requested quartet, full-energy agreement, decomposition closure and stop conditions remain governing.
