# Independent terminal review:16602

21 October 2026 UTC. **The attempt failed in native initial-guess construction, before molecular SCF. It supplies no energy, gradient or force-comparison result.** The decisive error is unavailable La atomic orbitals for the PAtom guess, not an observed memory failure or failed electronic convergence. This review used source/output text, saved fields and hashes only; no scientific replay or new job was run.

## Exact terminal evidence

Workspace: `workspaces/lanm_coupled_mechanics_20261017/lbt_omol_target_reference_v1/run/execution`.

| Artifact | SHA256 |
| --- | --- |
| RESULT.json | `881852372e543e2f2c721e1b162997cb348c2bc3ddbbbb7d6e910560f7b5c55f` |
| central_La/endpoint.out | `770979dbf59e185ea50ead6687278985f3bc0c4333e0c81604a0ed4756f990c0` |
| central_La/endpoint.out.execution.json | `d24e2e1c2134ead4a2b628eca2415bf8cff06cddcd858bf1a0f79077f4f4df8a` |
| NATIVE_EVIDENCE.json | `cfbe6e3a407d80e8fd54b86444c0d979869ca7d88b3fc8d0f5681e4afd2c53ab` |
| FORCE_COMPARISON.json | `28ecbe443ca13e944b62f6bf5e0fc9b21c2dc8a801f69cc48c57f177fc91c70c` |

Selected source/configuration/manifest, input/XYZ/external-basis and receipt artifact hashes remain exact. TARGET is `aa62db4f…`, CONFIG_RUN `1c9f22df…`, manifest `1b543ed7…`, input `d482826a…`, external basis `7389c56f…`. The final native output is161600bytes; scheduler stderr is empty because native stderr was captured in endpoint.out.

The result is `invalid`, with one native runner invocation/one attempted DFT call and zero other-engine, displaced or optimizer calls. The comparison is explicitly unavailable. The receipt records no launch error but `normal_termination=false`, `scf_converged=false`, `expected_engrad_complete=false`; no engrad exists. Its top-level ORCA return code is0 despite the internal MPI/GUESS abort. The existing completion checks correctly reject that misleading return code. The saved terminal event reports scheduler FAILED1:0,539s and24 allocated CPUs.

## Failure and method evidence

- Native line2154 and subsequent ranks report `orca_guess/guess_eht.cpp`, line712: **SCF orbitals not available for element57(La)**. Line2225 identifies `INITIAL GUESS: POLARIZED ATOMIC DENSITIES`. The ending explicitly says error termination in GUESS; MPI reports signal6/Aborted. These diagnostics identify the unavailable atomic-guess support in this execution. They do not distinguish a missing installation data asset from an unsupported entry in that backend's atomic library.
- No molecular SCF iteration or terminal residual block is present. RESULT's SCF-evidence exception for missing/ambiguous “Energy change” is therefore expected from this early failure, not evidence of hidden converged SCF or an admission-parser failure requiring repeated chemistry. There is no final electronic energy or analytic gradient to recover from the saved files.
- Startup completed in497.9s and printed maximum startup memory2522.2MB. Earlier matrix checks said the data fit. Allocation remained24MPI/96000MiB/MaxCore3100MB per rank. No out-of-memory or scheduler-kill diagnostic appears in the inspected output. These facts do not establish peak aggregate memory, but they give no basis for treating additional RAM as the remedy to the explicit La guess error.
- Native method settings reached ωB97M-V/libXC, UHF, charge−2, multiplicity1,1030 electrons,6239 orbital and8062 AuxJ functions, DIIS on/SOSCF off, MaxIter300 and the published integral thresholds. The earlier automatic SCNL and AutoTRAH-off warnings remain authoritative; ECP was labeled `def2-tzvpd.bas`, La46 at atom255. Successful execution of the beta rotation and a final electronic solution were not established.

## Recoverability and frozen controls

This is a native initialization compatibility failure, not evidence that the model/reference forces agree or disagree, nor that the requested electronic state cannot converge. An unchanged retry, larger MaxIter, looser tolerances or added memory is not supported by the observed cause.

A separately declared recovery could address the unavailable atomic guess while preserving this failed attempt. A documented alternative such as PModel would be an explicit initial-guess deviation, not exact reproduction of the failed PAtom deck or proof of the same SCF solution. Prior La PModel success establishes useful compatibility evidence, not success for this diffuse UKS calculation. Do not reuse the failed GBW as an admitted electronic seed or silently drop the beta rotation.

Keep the exact256 geometry/order/protonation, charge/multiplicity/electron count, external orbital basis and ECP, def2/J, ωB97M-V, executed self-consistent NL treatment, grids, DIIS/NOSOSCF/Tight/300/thresholds, fresh-directory policy and native admission requirements frozen except for a specifically declared supported guess change. Preserve the version limitation and unqualified physical state. Any future SP pair still requires its own consistent-NL declaration and actual central admission/cost review. No retry, setting change or additional call is authorized by this review itself.
