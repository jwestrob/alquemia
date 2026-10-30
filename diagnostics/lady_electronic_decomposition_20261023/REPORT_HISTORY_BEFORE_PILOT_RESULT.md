# Current result status: pilot unfinished after initialization OOM

21618 exhausted its96,000MiB allocation before the requested electronic states
were obtained. **No CT result is available.** Root submitted one memory-only
recovery21623 at256GiB/24CPUs onstandard_expand. Its7073s native cap retains the
original cumulative pilot budget; all electronic inputs and gates are unchanged.
The remaining three systems still require pilot admission.

Current evidence: `CURRENT.md`, `cp2k_memory_recovery_v1/RECOVERY_SPEC.md`,
`FAILED_ATTEMPT.json`, `RESOURCE_EVIDENCE.json`, actual workspace receipts.
Older status sections below are retained history.

---

# Bounded electronic diagnostic — pilot in progress

The active goal is the exact four-geometry La/Dy metal–whole-shell ALMO
comparison. No CT differential is available yet. **Compatibility pilot21618
has been submitted on24CPUs with standard_expand and a7200-second native cap.**
The remaining three are prepared and await actual pilot admission; they will
run concurrently only if the declared checks pass.

See `CURRENT.md`, `cp2k_inputs_v1/RUN_SPEC.md`, the immutable CONFIG and actual
workspace SUBMISSION/RUN_RESULT. The installation report below is historical;
its earlier paused-goal wording is superseded by the active goal continuation.

---

# La/Dy electronic decomposition: software compatibility status

Updated 23 October 2026. **The isolated CP2K2026.2 build is installed and passes startup/library checks.** No new electronic calculation has run. The CT differential is not yet known.

## Inputs ready

The four requested geometries, original native input decks and metal basis/ECP assets are copied without modification under `inputs/`. `INPUT_AUDIT.json` and `inputs/PROVENANCE.json` record the checks and original result locations. Nonmetal coordinates agree exactly within each La/Dy pair. The aquo systems have charge +3 and 88 explicit electrons; the acetate/water systems have charge +2 and 110. Metal and complete-shell fragments have effective singlet states. Physical Dy spin metadata must not replace the archived effective state.

## What the software investigation established

| Implementation | Finding for this exact task |
|---|---|
| Q-Chem | Documented ALMO and custom ECP support. No accessible licensed installation located; Jacob confirms he does not own a license. |
| Stock CP2K | Native ALMO and semilocal ECP integrals exist, but the stable 2026.2 ECP reader rejects the exact Dy core count of 55 with `Unknown Core State`. Inline potentials and later `ELEC_CONF` do not bypass that branch. Installing the stock binary would not fix it. |
| GPU4PySCF | Native ALMO-EDA2 exists, including PBE0 charged-fragment tests. Its EDA energy path explicitly asserts that no ECP is present. Generic GPU4PySCF ECP support does not make this EDA implementation compatible. |
| Qbics | Promising documented TSO/BLW implementation: constrained and full states, polarization/CT decomposition, PBE0, explicit Gaussian ECPs and a frozen-f Ce ECP47 aquo example. Exact Dy55 compatibility and state conventions remain unverified. Download requires registration; public license terms have not been established. |
| GAMESS-BLW | A documented BLW-DFT implementation exists in an author-supplied/in-house GAMESS version. A stock GAMESS installation alone has not been established as the required implementation. |

CP2K's initial charged-fragment concern is narrower than an outright ALMO prohibition. Its molecular initialization provides raw isolated-fragment and antisymmetrized frozen energies with nonperiodic boundaries. That supports a reconstruction route independently of the documented charged-sublattice warning. The separate Dy55 input-reader restriction remains decisive for the unmodified release. See `CP2K_PILOT_REVIEW.md` and `CP2K_COMPATIBILITY.md` for source evidence.

## Current route: isolated CP2K reader repair

Two independent source reviews identified a narrow repair: add Dy's existing 55-electron core occupation case to the ECP input reader. CP2K's downstream pseudopotential-state code already handles this core. The patch adds seven lines, preserves all ECP coefficients and the effective charge, and leaves the ALMO and integral equations unchanged. It does not construct a substitute electronic-structure method. See `cp2k_fix/INDEPENDENT_REVIEW.md` and `PATCH_SCOPE_REVIEW.md`.

The official CP2K 2026.2 source archive has been downloaded and its SHA256 verified against the publisher's release digest. The patch is applied only to the isolated workspace source; pristine source remains in the original archive. Actual before/after hashes match both reviews. `workspaces/lady_electronic_decomposition_20261023/cp2k_build_v1/PATCH_APPLIED.json` records the change.

The first attempt, job21396, downloaded the immutable official MPICH container successfully, then failed before compilation because the runtime image lacks CMake and pkg-config. It used1,197 seconds on24 allocated CPUs; zero molecular calls. The actual container contains the required quantum libraries, compiler, MPI and Ninja. A contained recovery supplies portable CMake3.28.3 and a separately pinned pkgconf build tool, reusing the exact saved source, image and quantum dependencies. Its driver is `cp2k_fix/BUILD_RECOVERY.py`, with actual receipts under `workspaces/lady_electronic_decomposition_20261023/cp2k_build_v2/` when present. The recovery retains24 CPUs under `standard_expand`, a two-hour process cap and version/loader checks. This is software setup, with no molecular evaluations.

Recovery job21440 completed in700 seconds on24 allocated CPUs (16,800 allocated CPU-seconds). Version and dynamic-library checks passed, and root independently verified startup without the build tools mounted. The binary SHA256 is `19e4c45c2920299863d0d74543b55ded033ab1e145cbba1cd5cfc3d660e23b32`. Both setup attempts together cost45,528 allocated CPU-seconds. See `INSTALLATION.md` for exact pins and invocation.

Build success does not establish scientific compatibility. Exact input translation, electron/state checks, the archived full-energy comparison and ALMO reconstruction remain required before accepting the La-aquo pilot or releasing the remaining three systems. The entire shell remains one fragment; native full SCF is required after constrained polarization. D4 and reference terms must be carried consistently. No result is yet available for the requested charge-transfer differential.

## Qbics fallback checks, if needed

Qbics calls its implementation TSO/BLW, not Q-Chem ALMO-EDA2. Its documented counterfactual is relevant, but naming similarity is not a compatibility test. Inspect the native source and license obtained through the official download route, then verify:

- The declared two fragments are M3+ and the entire ligand shell, with correct charge and electron-count interpretation.
- The frozen reference is the declared antisymmetrized fragment state and the constrained state suppresses metal–shell delocalization while allowing intrafragment response.
- PBE0, exact original ECP46/ECP55 coefficients, radial powers, orbital basis and spherical conventions are preserved.
- Its required unrestricted execution retains the archived effective singlet solution; it must not silently introduce a different broken-symmetry reference.
- No D3 correction is enabled. Its documented D3 option is not archived D4. Carry the original geometry-dependent D4 contributions separately and consistently in constrained/full totals.
- Raw and counterpoise-corrected terms are kept distinct and reconstruct the same declared total.

Only after compatibility is established does the authorized La-aquo pilot become executable. The other three systems remain conditional on that pilot and would run concurrently under `standard_expand`. No model fitting, new geometries or changed electronic protocol follow.

## Existing result, not a new decomposition

The archived monodentate water-to-acetate metal double difference is +0.5750525133 kcal/mol for frozen-f PBE0-D4 and -6.2358150700 kcal/mol for transferred AMOEBA (model minus reference: -6.8108675833 kcal/mol). These are conditional fixed-geometry energy contrasts, not binding free energies. No CT-suppressed, polarization or CT value is available from the current work.

## Recommendation

When the paused scientific goal resumes, verify the exact archived protocol through the bounded La-aquo compatibility pilot using this installed build. No user installation, registration or purchase is needed for this route. The Qbics registration question is superseded. If the repaired implementation fails the scientific checks, report the specific failure; do not change the Hamiltonian or automatically expand the experiment.

## Primary sources

- [CP2K stable ECP reader](https://github.com/cp2k/cp2k/blob/v2026.2/src/atom_types.F)
- [CP2K ALMO controls](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/ALMO_SCF.html)
- [GPU4PySCF native EDA implementation](https://github.com/pyscf/gpu4pyscf/blob/master/gpu4pyscf/properties/eda.py)
- [Qbics constrained-state documentation](https://qbics.info/doc/keywords/scfguess.html)
- [Qbics decomposition definitions](https://qbics.info/doc/keywords/eda.html)
- [Qbics pseudopotentials](https://qbics.info/doc/keywords/pseudopotential.html)
- [Qbics official download and registration](https://www.qbics.info/home/qbics-download.html)

No registration, external contact, molecular calculation or job cancellation has been performed. Software setup is complete in the isolated workspace; existing scientific assets and unrelated jobs are preserved. The goal remains paused, as confirmed after installation.
