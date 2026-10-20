# Independent endpoint review — 15321

**Verdict: admitted under the unchanged native contract; no blocker to the saved collector.** The one requested `new_La` endpoint completed normally with a full 256-atom analytic gradient. This is native availability, not numerical force qualification, all-criteria SCF convergence or proof of a stable electronic state.

## Exact evidence

Workspace: `workspaces/lanm_coupled_mechanics_20261017/lbt_7ccn_incremental_reference_v1/run/execution/`.

| Artifact | SHA-256 |
| --- | --- |
| `RESULT.json` | `3a173d81072ab269922d17db1596c130e8132bd2f803614b00d2fe59455c10de` |
| Frozen source | `bf75f502a2e3c53e51329cfd59b9ff34f661d85acc0ed24483f51868bbf194cd` |
| Run configuration | `155d0d6ebc0b665a6ce287e519092d1b25a1ed55f5d4eba779189502fe45a996` |
| Native receipt | `6583718c6fbe8aa5aa1843b2ab1fefc9e877977c9990f484560f00a04be1b9d9` |
| Native output | `e514f1c3f9b2f2d20a67b43b7ae234e53fb6bbb1e2ee176c38fd512e7522291b` |
| Native gradient | `8ee8c82dc21f1df48ca38d8ecf8bb0d57bd91fcf24b4172a80f0797eab0702a2` |
| Run manifest | `d96eda727e33962a1195d2ebaa7862587ae71ac6641819104b2b584288044d7d` |

## Admission and scope

- The run binds actual preflight15319 `e74495a8` and execution decision `05b4ad31`. Its configuration differs from the frozen preflight configuration only in those two completed bindings. Scientific binding, input evidence, both solvent admissions and the old native row are exact saved matches.
- The sole task matches the actual preflight signature: exact terminal POINT `eb96df29`, physical inventory and coordinates, unchanged native template, charge −2, La physical/effective singlet, ECP46 and 1030 explicit electrons. Original MaxIter125 and vacuum/PModel recipe remain. The printed atomic-guess multiplicity 2 does not change the molecular singlet declaration.
- Exactly one native runner invocation and one DFT endpoint are recorded. All new MACE, ALPB, GFN, MM, optimization, geometry and other-metal counters remain zero. The old reference14880 and both solvents are reused.
- ORCA6.1.1 returned zero, reported SCF convergence and terminated normally. The receipt marks the expected gradient complete. The saved parsed gradient is 256×3; the native gradient section contains all 768 components. Final single-point and engrad energy are printed as −6765.706948267616 Eh.
- Runtime used 12 MPI ranks, one thread per rank, one worker, 48000 MiB allocated memory, 25% reserve and MaxCore3100. Native receipt duration is 24858.128134 seconds, finishing at 2026-10-20T09:07:22.491989+00:00. Scheduler accounting remains authoritative for allocated cost.

## Printed SCF evidence

The new endpoint explicitly reports **Energy Check signals convergence**, followed by convergence after **15 cycles**. `ConvCheckMode=Total+1el-Energy` and `ConvForced=0` are unchanged.

| Printed terminal criterion | Value | Printed tolerance | Pass |
| --- | ---: | ---: | --- |
| Energy change | −5.5297e−10 | 1e−9 | Yes |
| Maximum density change | 5.8322e−5 | 1e−8 | No |
| RMS density change | 3.7163e−7 | 1e−9 | No |
| Orbital gradient | 5.0204e−6 | 2e−6 | No |
| Orbital rotation | 8.2374e−6 | 2e−6 | No |
| DIIS error | 1.1692e−3 | 1e−8 | No |

Old14880 had the opposite energy/orbital-gradient flag pattern: energy failed its printed tolerance while orbital gradient passed. Its exact original evidence is retained. Neither endpoint passes all printed criteria. The new small-gap warning (−0.019) remains explicit. Final gridx3 integration completed; its printed scalar correction is not a force-error estimate.

## Interpretation limits

Both rows are available, but this run's comparison remains unavailable by design pending the saved collector. No work difference was recalculated here. Full force reports retain 256 physical IDs and descriptive comparisons; `physical_state_qualified` and `force_consistency_qualified` remain false. Normal termination does not resolve vacuum-dianion binding, frozen-f/ECP, conditional terminal chemistry, post-added fixed-solvent or native-gradient accuracy limitations. No retry, new endpoint or optimizer extension follows from this admission.

Review verified 467 explicit file/hash pairs and literal saved bindings, inspected native text and gradient dimensions, and performed no scientific replay, scientific imports, molecular calls, submissions or source changes.
