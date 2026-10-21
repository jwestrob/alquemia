# Independent startup review:16602

Observed21October2026 at03:12:39UTC. **The native process has passed initial input parsing, with the declared24-rank allocation and exact prepared inputs. A material method warning requires root review: ORCA automatically switches the EnGrad calculation to self-consistent SCNL.** No job, source or configuration was changed by this reviewer.

## Actual allocation and identity

Saved SCONTROL reports RUNNING on `node-64-768g-13`, standard/standard_expand,24 tasks/24CPUs,96000MiB, one node and no restart. The runtime sidecar records one worker,24 local MPI slots,25% memory reserve and MaxCore3100 decimalMB per rank. The runtime input preserves the prepared scientific deck and adds only the allocation-derived MaxCore and24-process PAL block.

The run task signature exactly matches actual preparation16601, including input, XYZ and adjacent external-basis hashes. Run configuration differs from preparation only in the two released future pins. Source, decision, configuration, manifest, input evidence and task-file hashes were checked. Input evidence remains byte-identical to preparation.

| Artifact | SHA256 |
| --- | --- |
| CONFIG_RUN.json | `1c9f22df1c9ac096568dccac7689a5e6c9f8af4d8f8ae2bd0a3d589aa9ca0c83` |
| ROOT_EXECUTION_DECISION.json | `7c7e0dcc94d2ff1d862113f362d4ea21b265eafeb3846ec5a03e5f1447bcde5f` |
| execution/START.json | `e772eae8c3d8b36c7d1f9f684f7b29af5f9134e5947c711022b15bad2113dee1` |
| execution/MANIFEST.json | `1b543ed7c470cd9a8df7b0a7cfa3329cf6c8f30be862f896ec15049c9a44d43e` |
| execution/SCONTROL_JOB.txt | `05a84f2c8d4fbef272f42164eb7e5517d015e941c009ef6a71b54118bd6695c1` |
| central_La/endpoint.runtime.inp | `d028084243c29e768fa5ab609614de7eb4ad3504ba0552800f51e55a93becd23` |
| central_La/endpoint.runtime.json | `f7daac838334126e20573a6df5df16bd680785ae6ad30a40cf424122fc14c134` |

## Executed method evidence and correction

The changing `execution/central_La/endpoint.out` was103742bytes, last modified03:12:02.433400UTC, with1929lines at observation. No final immutable output hash is asserted.

- Line54 identifies ORCA6.1.1. Lines203–210 identify the orbital basis read from `def2-tzvpd.bas` and auxiliary def2/J. Line1835 identifies **La ECP `def2-tzvpd.bas`,46 replaced core electrons,lmax3**, with atom255 assigned at1837. This supersedes the uncertain source-only expectation that the executed ECP would necessarily carry a built-in label. The prepared file/deck is unchanged; detailed coefficient review remains pending.
- Lines218–220 state that the libXC B97M-V/wB97M-V variant is activated and manually given DFX/HFX/RangeSep scalings are ignored. The exact prepared deck supplies no such manual scaling overrides.
- **Lines223–227 explicitly state that analytic DFT-NL gradients require the self-consistent method, switch to SCNL, and switch AutoTRAH off.** Therefore the declared/input-evidence string “post-SCF VV10” is not an accurate description of the actual EnGrad treatment. This is observed native behavior under the published deck, not an input substitution or parser failure. It does not establish whether historical ORCA6.0.0 behaved identically.
- The latest tail reports that four symmetric matrices of dimension6239 fit within the current MaxCore allowance. No SCF iteration, final state header, executed PAtom/rotation evidence, converged result or native receipt was available at this observation. UKS/PAtom/rotation are present in the exact input; execution of those stages is still pending.

Completion watcher42618 and health watcher42619 have receipts for this exact job/manifest. The health receipt at03:12:23UTC is `watching`, with no explicit failure or alert and no terminal endpoint. These receipts establish registered observation, not useful numerical progress or convergence.

## Decision limits

Root was notified immediately of SCNL and the actual ECP label. There is no observed parsing/resource stop trigger. **Continue the exact published EnGrad input unchanged:** its empirical6.1.1 energy/gradient remains informative about what that recipe actually executes at the fixed consumed geometry. Stopping solely because the earlier default interpretation was wrong would discard that bounded question without demonstrating that the published input is wrong. This recommendation does not authorize altered settings or additional calls.

The result must acknowledge the executed SCNL treatment and cannot be admitted as the originally described post-SCF calculation. Preserve the original source/SPEC expectation and this correction separately. No displaced call is released. Any future singlepoint comparison must explicitly resolve EnGrad-versus-SP nonlocal treatment before launch, rather than silently comparing different Hamiltonians. Full method/state/force admission and measured cost remain pending; the failed PBE0 gate is unchanged.
