# ExaA donor-response prelaunch review

Independent static source review; no molecular imports, tests, submissions or calculations.

## Scientific scope and metadata binding

No blocking mismatch found in `PLAN.md`, `reference/PREPARE.py` or `reference/RUN.sbatch` against the actual frozen inventory `757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343`.

- The plan declares sample3 target/sample0 control, exact 192-atom origins and saved adaptive_La geometries, both metals: eight vacuum endpoints. The common path is explicitly distinguished from the classifier's separate selected geometries.
- PREPARE uses the actual inventory field names, copies exact pinned XYZ bytes and state/admission evidence, and binds the original native runtime/executable. It does not introduce geometry changes, a solvent term, new electronic Hamiltonian, fitted boundary or automatic restart.
- Inputs retain r2SCAN-3c/NoAutostart/DefGrid3/TightSCF with recorded Ca−3/La−2 singlets. Electronic counts and geometry admission are left to the declared real Slurm preflight rather than guessed from atom count.
- The wrapper requests one standard node and 12 MPI slots with numerical libraries at one thread per rank. It specifies neither a memory/time override nor a priority/account/QOS workaround.
- The first two scouts are the target La origin/proposal, each in its own 12-rank job. Their first result is an individual La work; it does not by itself establish the paired Ca-minus-La response.
- The plan preserves all eight cells, real failures and unknown electronic-state stability. Its 1kcal/mol interpretation scale is identified as a declared resolution floor, not a proven reference uncertainty or fitted scientific acceptance gate.

## Execution adapter

Completed focused source review of `reference/REFERENCE.py`. **No launch-blocking issue found.** Actual source JSON and a real completed native receipt were read to check helper/schema compatibility; no scientific functions were imported or executed.

- Source/admission/map/state/proposal joins match the actual inventory. The preflight derives electronic counts including sulfur, verifies singlet parity, exact saved XYZ bytes and paired coordinates, and requires eight separate output directories. It checks the actual allocation before native work.
- The corrected TRAH collector receives the pinned `OLD_REFERENCE.py` module directly; it does not accidentally reload the new adapter through the old loader. Its real missing-summary water fixture and existing population/parser checks execute only on Slurm.
- Sorted nonblocking locks cover only each requested endpoint, with execution and collection inside the lock. The two distinct scout tasks therefore need not serialize. The unchanged native runner writes only task-local files and refuses partial/unverified outputs or stale runtime files; no automatic retry or overwritten scientific result is introduced.
- Each attempt directory contains the actual job/action/task identity. Collection reads only the selected endpoints and explicitly records the full eight-cell declaration. A failed required endpoint remains failed/missing; no energy substitution or hidden successful fallback occurs. Retained GBW/populations and actual ending-solver evidence follow the existing repaired collector.
- The new receipt checks match actual native artifact keys and require executable, runner, renderer, template, coordinates and 12×1 parallelism. Collection remains a conditional diagnostic: convergence is not electronic-state stability, response accuracy or classifier validation.

Reviewed file hashes:

| File | SHA256 |
|---|---|
| PLAN.md | `87c1bfe98e1a496ed6a8da6ec9558acd85c7cd946889b34ea46f8eccabcea520` |
| reference/PREPARE.py | `40f1df507d37b5cf3222468086533c4237c04521e4c17f89f189940d89059dcb` |
| reference/REFERENCE.py | `7ec03787224f05785741c611a3700e8ad470e0e7c7f693dc876dd6712238e894` |
| reference/RUN.sbatch | `df688249568168254573fe875596f986b61650453991a13b2e29b8abb19d081d` |

This is static review, not an executed runtime/preflight pass. Root owns actual submissions and acceptance. No implementation files were changed by this reviewer.
