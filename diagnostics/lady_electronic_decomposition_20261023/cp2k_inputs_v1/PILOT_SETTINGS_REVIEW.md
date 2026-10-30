# Exact-release ALMO pilot settings review

2026-10-23. Independent source/documentation review only. No quantum execution, parser execution, installation, submission or numerical convergence experiment was performed. Root constructs and binds the final deck.

**Finding:** the inspected release supports the proposed native two-fragment molecular-guess → frozen → block-polarized → FULL_SCF path. No hard combined GAPW/PBE0/semilocal-ECP incompatibility was found in the inspected call path. This is source support, not a tested combination or a promise of convergence/0.1 kcal/mol agreement. The Dy55 reader patch is separately reviewed; source/input/build identity must remain explicit.

## Authority

Actual source root: `workspaces/lady_electronic_decomposition_20261023/cp2k_build_v1/cp2k-2026.2`.

The downloaded v2026.2 release archive has SHA256 `f9bd86f580f57a53a0768c0045d1417f9f9a1d66d851ed7f662f496200043373`; its `REVISION` says `git:c92cc08`. These actual archive/files govern this review, rather than assuming byte identity with an earlier master/tag citation. The official versioned manual corroborates interface meanings.

| Actual source file | SHA256 |
|---|---|
| `src/almo_scf.F` | `24d543aa07f902d1876153a7d473e4e980306fd1c44574e03aab51652944f426` |
| `src/input_cp2k_almo.F` | `d661bc2549e6bb45cddedf51f0de197aa727f54fb1f73008e74c94a039643fab` |
| `src/almo_scf_optimizer.F` | `01ce8e4ca6bd9f61cda35a061edab8a94600959d27221445f952bc519b3f8dc3` |
| `src/mscfg_methods.F` | `03480b825354ec80e61c4a85ee4b7db6faa0d99c4f8eb99dc36c18001a726d37` |
| `src/qs_ks_methods.F` | `87c2f9e7055e467695f40d92a389d8252e2d940d81daa0f4878ee7778d8791fa` |
| `src/cp_control_utils.F` | `da04a34bb800c6ccc9314a9e434abca7b0f1482a250f1672bb0cd27adff6aa39` |
| `src/topology_coordinate_util.F` | `c543a9a21539f1834f667ab60cdd4f4f4e951f2c646a82e18ec26f2d4a2f3807` |
| `src/pw/mt_util.F` | `33b9978bd6bafda8b03a32806ac2117e0bb0597f983efcecf59e769368774b38` |

## Geometry, groups and isolated boundaries

Keep the literal La-aquo coordinates and ordering: metal atom 1; entire water shell atoms 2–25. Append only CP2K's molecule name/index columns in the native COORD block:

```text
&COORD
  La 0.000000000000 0.000000000000 0.000000000000 METAL 1
  # Every remaining original O/H coordinate row ends in: SHELL 2
&END COORD
&TOPOLOGY
  CONNECTIVITY USER
  &CENTER_COORDINATES F
  &END CENTER_COORDINATES
&END TOPOLOGY
&CELL
  ABC [angstrom] 30.0 30.0 30.0
  PERIODIC NONE
&END CELL
```

This is two groups, not nine molecules. `regtest-almo-1/almo-singlet-charge.inp` demonstrates USER grouping and explicit charged FRAGMENT sections. No fictitious bonds are introduced. Verify the actual output reports exactly two domains and occupancies 4 and 40 doubly occupied orbitals (88 total electrons), charges +3 and 0, and all-neighbor/full delocalization.

**Negative coordinates need no silent translation.** In actual `topology_coordinate_util.F:315–369`, the nonperiodic check compares coordinate extents with cell lengths; with centering false, `dims=0` and coordinates are copied unchanged. Actual `tests/QS/regtest-hybrid-1/H2O-hybrid-pbe0.inp` already uses negative coordinates with GAPW/PBE0/MT and no CENTER_COORDINATES. MT constructs a relative-coordinate nearest-image Green function (`mt_util.F:150` onward); placing the cluster across the grid origin is not by itself an input-domain failure. No alternative cell-origin keyword is needed here.

Use `DFT/POISSON/PERIODIC NONE` and `POISSON_SOLVER MT` as well. The 30 Å cube is one deliberately generous prospective choice for these compact clusters, not a demonstrated density-tail bound. MT requires a cell larger than twice the relevant cluster/density extent. Do not substitute WAVELET's face-density criterion or a periodic charged-cell background. An actual cell/accuracy discrepancy remains a failed pilot, not permission for an automatic cell sweep. [2026.2 Poisson documentation](https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/POISSON.html), [optional centering](https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/TOPOLOGY/CENTER_COORDINATES.html)

## One fixed numerical realization

Propose these values prospectively, unchanged across the quartet if the pilot is admitted. They are numerical choices, not measured error bars:

```text
# Under DFT
CHARGE 3
MULTIPLICITY 1
UKS F
&MGRID
  CUTOFF [Ry] 1200
  REL_CUTOFF [Ry] 100
  NGRIDS 5
&END MGRID
&POISSON
  PERIODIC NONE
  POISSON_SOLVER MT
&END POISSON
&QS
  METHOD GAPW
  ALMO_SCF T
  EPS_DEFAULT 1.0E-12
  EPS_PGF_ORB 1.0E-10
  EPSFIT 1.0E-6
  EPSSVD 1.0E-10
  EPSISO 1.0E-14
  EPSRHO0 1.0E-10
&END QS
&SCF
  SCF_GUESS ATOMIC
  EPS_SCF 1.0E-10
  MAX_SCF 300
  IGNORE_CONVERGENCE_FAILURE F
  &OT
    MINIMIZER CG
    PRECONDITIONER FULL_ALL
  &END OT
&END SCF
```

Add `RADIAL_GRID 200` and `LEBEDEV_GRID 590` to each KIND. 590 is an actual supported Lebedev grid size. Preserve the exact archived orbital basis and semilocal ECP under `ECP_SEMI_LOCAL T`; all-electron O/H use `POTENTIAL ALL`. The later mixed shell also requires exact C functions. Do not import an auxiliary RI/ADMM basis as an orbital basis.

The staged conversion uses `BASIS_SET LADY-ARCHIVED` for every kind, `POTENTIAL ECP ECP-LADY-La46` for La (`ECP-LADY-Dy55` for Dy), and `POTENTIAL ALL` for H/C/O. Its separate `assets/PROVENANCE.json` and `CONVERSION.md` govern coefficient mapping. The archived metal coefficients are explicit; original light-atom library coefficients were not printed, so H/C/O named-library/BSE equivalence must remain labeled as that evidence level rather than literal archived-coefficient verification. This review has not rerun the conversion.

`EPS_DEFAULT` alone is insufficient to describe Gaussian screening: `cp_control_utils.F:1009–1028` sets `EPS_PGF_ORB = sqrt(EPS_DEFAULT)` unless explicitly overridden. The proposed override avoids silently accepting 1e-6 orbital screening. GAPW's hard/soft split has independent `EPSFIT/EPSSVD/EPSISO/EPSRHO0` controls; making EPSFIT smaller also makes the grid burden harder. The high cutoff and denser atom-centered grids are intended to address both parts. No parameter guarantees the requested cross-program agreement. The 30 Å/1200 Ry grid has material memory cost; do not reduce it silently to fit an allocation. Bind sufficient scheduler memory and actual usage.

The ordinary SCF block controls the **two isolated fragment calculations**. It is not the final ALMO stopping control. `mscfg_methods.F:174–319` retains the same DFT/FORCE_EVAL and cell, turns ALMO off internally, passes fragment charges/multiplicities and calls `qs_energies`. The comments about future fast/dirty shrinking are not active changes to this path. No OUTER_SCF loop, smearing or restart guess is requested.

## ALMO controls and complete delocalization

```text
&ALMO_SCF
  ALMO_SCF_GUESS MOLECULAR
  ALMO_EXTRAPOLATION_ORDER 0
  XALMO_EXTRAPOLATION_ORDER 0
  ALMO_ALGORITHM DIAG
  XALMO_ALGORITHM PCG
  DELOCALIZE_METHOD FULL_SCF
  MO_OVERLAP_INV_ALG DENSE_CHOLESKY
  EPS_FILTER 1.0E-12
  CONSTRUCT_NLMOS F
  &ALMO_OPTIMIZER_DIIS
    EPS_ERROR 1.0E-9
    MAX_ITER 300
    N_DIIS 8
    EPS_ERROR_EARLY -1
    MAX_ITER_EARLY -1
  &END ALMO_OPTIMIZER_DIIS
  &XALMO_OPTIMIZER_PCG
    EPS_ERROR 1.0E-9
    MAX_ITER 300
    MAX_ITER_OUTER_LOOP 0
    CONJUGATOR HAGER_ZHANG
    PRECONDITIONER FULL
    LIN_SEARCH_EPS_ERROR 1.0E-3
    LIN_SEARCH_STEP_SIZE_GUESS 0.1
    EPS_ERROR_EARLY -1
    MAX_ITER_EARLY -1
  &END XALMO_OPTIMIZER_PCG
  &ANALYSIS F
  &END ANALYSIS
  &FRAGMENT 1
    CHARGE 3
    MULTIPLICITY 1
  &END FRAGMENT
  &FRAGMENT 2
    CHARGE 0
    MULTIPLICITY 1
  &END FRAGMENT
&END ALMO_SCF
```

The exact names are `MO_OVERLAP_INV_ALG`, `DENSE_CHOLESKY`, and `DELOCALIZE_METHOD FULL_SCF`. In `almo_scf.F:1624–1655`, FULL_SCF explicitly uses no spatial quench, `perturbation_only=.FALSE.`, and the fully delocalized special case. It rejects truncated virtual spaces; keep the default full virtual space. `XALMO_ALGORITHM DIAG` is advertised generally but is not accepted by this FULL_SCF branch, hence explicit PCG. Outer-loop zero avoids silently multiplying the iteration cap. Failure to reach the declared error threshold is fatal with early stopping disabled. [2026.2 ALMO interface](https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/ALMO_SCF.html)

## Exact PBE0 exchange path

```text
&XC
  &XC_FUNCTIONAL NO_SHORTCUT
    &PBE
      SCALE_X 0.75
      SCALE_C 1.0
    &END PBE
  &END XC_FUNCTIONAL
  &HF
    FRACTION 0.25
    &INTERACTION_POTENTIAL
      POTENTIAL_TYPE COULOMB
    &END INTERACTION_POTENTIAL
    &SCREENING
      EPS_SCHWARZ 1.0E-12
      SCREEN_ON_INITIAL_P F
    &END SCREENING
    &MEMORY
      MAX_MEMORY 512
    &END MEMORY
  &END HF
&END XC
```

The 512 MB HFX cache is per MPI process; allocation must include its aggregate and the grids. This is a storage choice, not a changed interaction. No ADMM, ACE, finite-range Coulomb, RI-HFX, COSX surrogate or range-separated functional is proposed. CP2K and archived ORCA RIJCOSX are different numerical realizations of the same intended PBE0, so full-energy agreement remains an actual test.

The release regression `H2O-hybrid-pbe0.inp` combines GAPW, MT and these PBE scaling/HF fractions. Semilocal-ECP GAPW has a separate native regression. ALMO sends its density to the common KS evaluator (`almo_scf_qs.F`, `almo_dm_to_almo_ks`), which includes GAPW local-density terms and ordinary full HFX (`qs_ks_methods.F:778–812`); no ALMO-specific hybrid exclusion was found. These source paths do not replace runtime qualification.

## Raw output, frozen term and admission

Use `GLOBAL/PRINT_LEVEL MEDIUM` and `RUN_TYPE ENERGY`; LOW also retains the ALMO lines below because they write directly to the active logger. Keep full stdout/stderr, not only a summary. A fresh output directory and atomic fragment SCF guess prevent reading unrelated wavefunctions.

| Quantity | Actual native source/output |
|---|---|
| Each isolated fragment energy | Its ordinary converged QS/SCF result in the fragment section |
| `R = E_metal + E_shell` | `Single-molecule energy:` in `almo_scf.F:792–795` |
| Frozen absolute `F` | `Energy of the initial guess:` immediately after native antisymmetrized molecular density evaluation |
| Polarized absolute `P` | Last converged `ALMO SCF DIIS` row; error norm checked by native DIIS exit condition |
| Full absolute `U` | Converged full-PCG result and final `ENERGY| Total FORCE_EVAL ( QS ) energy` |

**Parser correction:** `ENERGY OF BLOCK-DIAGONAL ALMOs`, `ENERGY LOWERING` and `CORRECTED ENERGY` are emitted by the perturbative reporting path. They must not be required as the sole source of P/U for FULL_SCF. `wrap_up_xalmo_scf`'s nonperturbative branch updates the final energy directly. The chosen block DIIS trace supplies P with ten decimal places.

Native quantities support combined `FRZ=F-R`, `POL=P-F`, `CT=U-P`, with `FRZ+POL+CT=U-R`. I found no implemented native isolated electrostatic/Pauli split in this path. Report that split unavailable; do not reinterpret Hartree/XC components as an EDA partition or enable charged SUBLATTICE analysis to invent one. Whole-shell internal delocalization remains permitted.

Preserve all fragment/block/full convergence and occupation records. Full electronic energy versus archived ORCA electronic energy, and reconstructed PBE0-D4 total versus archived total, can use a prospectively fixed 0.1 kcal/mol agreement target if root adopts it. Reuse the exact archived per-geometry D4 once in the separate ledger; it cancels between F/P/U at fixed geometry. Do not enable a different empirical dispersion correction. Arithmetic component closure, ten-decimal printing, and small SCF residuals do not independently prove that every component has 0.1 kcal/mol accuracy.

Any native failure, unmatched effective state, basis/ECP discrepancy, incomplete raw-energy chain, or failed declared agreement remains a failed pilot. No automatic restart, alternative optimizer, grid sweep or modified decomposition follows from this settings review.
