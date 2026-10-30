# Latest checkpoint — recovery running charged-fragment SCF


Latest direct observation:21623 is live and has progressed into the actual
charged-La fragment SCF (8electrons/4occupied orbitals). First reported OT energy
and line-search step are present; this is not yet fragment convergence or an
ALMO result. PeakRSS about250,259,184KiB remains below262,144MiB grant. Exact
reported rows and live scheduler/usage snapshots are in the task's PROGRESS.json.
The initial failed21618 callback was actually received and ROOT_ACK.json saved.

The exact La-aquo pilot's initial attempt **21618 ended OUT_OF_MEMORY**,126s,
24CPUs,96,000MiB, before any requested fragment/ALMO energy completed. Raw output,
failed parser result, scheduler evidence and3,024allocatedCPU-s remain preserved.
This did not test the charge-transfer differential or establish a method defect.

Root reviewed and submitted **21623**, one memory-only recovery: shared
standard/standard_expand,24CPUs (12MPI×2OMP),**256GiB**. Exact electronic decks,
assets, binary, parser, thresholds and geometry are unchanged. Native cap7073s
conservatively deducts127s from the original7200s pilot budget. No further
resource escalation/restart follows automatically. See
`cp2k_memory_recovery_v1/RECOVERY_SPEC.md` and its independent reviews.

Current configSHA256
`7447297cead20630b1d2cb1aaa0e097a22948b5eb1a215afa91dc66b25ad2d29`.
Actual taskroot:
`workspaces/lady_electronic_decomposition_20261023/cp2k_quartet_memory_v1/`.
Pilot watcherPID639887, event `lady-cp2k-La-aquo-memory-recovery-v1`.
**Remaining three are unreleased** until actual pilot admission; then run them
concurrently with these exact prepared settings. No duplicate submission,
changed electronic protocol, fit, geometry generation or protein campaign.

The late21618 completion callback refers to the preserved failed allocation;
it does not authorize a duplicate retry. Read actual21623 state/receipts next.
`cp2k_inputs_v1/COLLECT.py` is prepared source-only for the eventual admitted
quartet; it has not produced a scientific result.

---

# Current software setup — 23 October 2026

Jacob asked to keep looking for a viable alternative and whether he needs to install it. **The isolated CP2K2026.2 installation is complete and passes startup/library checks.** No user installation, registration or purchase is needed for this route. The earlier Qbics registration question is superseded.

## Completed software setup

Slurm **21396 failed before compilation** after successfully downloading the pinned container. Actual image inspection found the expected dependency view and quantum libraries intact, but CMake and pkg-config absent from the runtime image. Root acknowledged the completion receipt in `ROOT_ACK.json`. Cost: 1,197 elapsed seconds, 28,728 allocated CPU-seconds; peak batch RSS9,607,432K. Zero molecular calls.

The original attempt, image and source are preserved in:

`workspaces/lady_electronic_decomposition_20261023/cp2k_build_v1/`

- `SUBMISSION.json`
- `BUILD_RESULT.json`
- `completion_event.json`
- `pull.log`, `configure_build.log`, `slurm-21396.out`

The recovery **completed as Slurm21440**, exit0:0,700 seconds on24 CPUs under `standard_expand` (16,800 allocated CPU-seconds; peak batch RSS5,713,196K). Its completion receipt is acknowledged in `workspaces/lady_electronic_decomposition_20261023/cp2k_build_v2/ROOT_ACK.json`. It reused the exact downloaded image and patched source, supplying portable CMake3.28.3 and pinned pkgconf build tools. Requested MPI/OpenMP, Libint, Libxc, FFTW, ScaLAPACK and D4 support are enabled. Native version and loader checks passed; root independently repeated the version check without build-tool bindings. All setup jobs are terminal; do not repeat them.

Executable: `workspaces/lady_electronic_decomposition_20261023/cp2k_build_v2/build/bin/cp2k.psmp`. SHA256: `19e4c45c2920299863d0d74543b55ded033ab1e145cbba1cd5cfc3d660e23b32`. It runs inside the retained pinned MPICH container. See `INSTALLATION.md` for the exact invocation and provenance. No production default was changed.

## Exact repair and remaining gate

The official CP2K2026.2 source archive digest is verified. The reviewed seven-line Dy55 default-occupation input-reader patch is applied only to its isolated workspace copy. Original archive and original reader file are retained. `PATCH_APPLIED.json` records exact source and patch hashes. ECP coefficients, core count, effective charge and ALMO equations are unchanged. Two independent source reviews passed.

The immutable official MPICH container supplies the runtime libraries. The separate executable passes version and dynamic-library checks; this does not establish scientific qualification. The independent source audit confirmed that the seven-line reader repair is the only change across all1,585 checked compiler-input files.

## Active goal: bounded quartet, pilot preparing

The latest goal continuation reactivated the explicit quartet objective. The
previous software-only paused status is historical and no longer a hold.
`cp2k_inputs_v1/` now contains all four exact native decks, converted orbital/ECP
assets, the saved reference energy ledger, and a fixed RUN_SPEC. Geometry tokens
and paired shells are unchanged. No molecular call has yet been launched.

Root is binding the one La-aquo compatibility pilot on24CPUs (12MPI×2OMP),
standard_expand, native process cap7200seconds. Per-cell cross-code energy gate
is1e-4Eh, with native convergence and decomposition checks separately required.
Only an admitted actual pilot releases the remaining three concurrently. An
incompatible pilot stops the bounded diagnostic; no automatic numerical retry,
protocol sweep, fitting, new geometry or protein campaign is authorized.

Read `cp2k_inputs_v1/RUN_SPEC.md`, actual CONFIG/SUBMISSION/RUN_RESULT and native
outputs before acting. Installation success is not scientific qualification.

### Pilot submitted

Root submitted **21618** for aquo8_central_La only,24CPUs on standard_expand,
12MPI×2OMP, native cap7200s. Exact config SHA256
`075aaf39b904651e0cd406fc127aa23523986c344b294fda1b61ebeb554db423`.
Actual SUBMISSION and watcherPID266404 are under
`workspaces/lady_electronic_decomposition_20261023/cp2k_quartet_v1/tasks/aquo8_central_La/`.
Completion event `lady-cp2k-La-aquo-compatibility-pilot-v1` is armed.
All four decks passed native syntax-only checks. This does not admit the science.
Do not resubmit, modify bound inputs, or release remaining three until root
reviews actual pilot output and writes an explicit PILOT_ADMISSION.json.

Actual startup reviewed:21618 RUNNING on node-64-768g-14, native12MPI banner,391 sphericalAOs and La46 ECP present. No energy/decomposition is admitted yet. STARTUP_REVIEW retains measured batch memory/CPU snapshot; allocation size is not utilization.

Recoverystartup actual:21623 RUNNING, grant262144MiB,12MPI; observedRSS108865416KiB already exceeds the previous96,000MiB cap. No requested ALMOresult admitted. Collector source reviewed and ledgerpath corrected for resourceconfiguration; finalSHA813b6e13a6afbd5d7137bface4e81bfff2bad257b1acb967f8d57c8e5ceb0366; unexecuted.
