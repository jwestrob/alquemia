# Independent static review: resource-only recovery

21 October 2026 UTC. **PASS for one-CPU zero-call preparation.** No static source or launcher blocker was found. Actual preparation and a separately bound root decision remain required before the single recovery endpoint. Review used text/diff, JSON and hashes only; no scientific imports, numerical tests or jobs.

## Frozen scope

| Artifact | SHA256 |
| --- | --- |
| TARGET.py | `5cad26238cc71e676a5af341ea8e0fe2a3ef4d7476f3a3d0503879c27ecc35c3` |
| SPEC.md | `29b3261699b11f8cbd905a6614765a9bc87db8205b094ab14e996f774e6033b6` |
| CONFIG_PREPARE.json | `32c8501b73d00df7affa6bf0a6d86937fe0ca2cf0b0028f36c5aa6bb37f96c9a` |
| preparation/RUN.sbatch | `99453665c44b014798f4d940fe3c63bc39f94aa9a30b6e349250543be744314a` |
| run/RUN.sbatch | `2ac03ea1617f64ed98981dede0bffb1a49742ed54bf893302b062189223e6d00` |

All19 direct configuration file pins match. Future preparation/result decision pins remain null. READY was not present at this observation; this review binds the actual inspected source/configuration and launchers above.

## Findings

- The new adapter fixes OOM16611's result `0fc5e6a1…`, source `27ac4e80…`, configuration `b654d846…`, input `4359537c…` and XYZ `55e07d81…`, while retaining the full failed16602 authority chain. Its guards check both exact failed-parent identities and the unchanged physical/electronic state.
- Generated scientific template, XYZ and external basis must match16611 by hash. The older exact PAtom→PModel relationship is also retained. No geometry, PModel/rotation, functional/SCNL, basis/ECP, solvent, grid, SCF threshold/limit or force-comparison arithmetic change appears in the diff. The only native runtime changes are the separately declared process count and allocation-derived memory.
- The run requires one node,12 tasks/MPI slots, two allocated CPUs per task,24 allocated CPUs and96000MiB. All five numerical-thread environment variables must equal1. The unchanged renderer must return one worker,12 ranks,24 allocated CPUs,12 slots,25% reserve and MaxCore6200MB. The terminal receipt admission also checks the actual rank/memory policy. CPUs-per-task is not substituted for MPI ranks or thread count.
- The inspected launchers implement exactly one-CPU preparation and12×2-CPU native allocation, with standard/standard_expand, no memory/time overrides, no requeue and explicit single-thread variables. Their explicit Python invocation accepts the script as its first argument. Preparation does not call the native allocation guard, so it does not require12 slots itself.
- Preparation/result/task-signature and root-decision gates precede the sole native invocation in a fresh directory. No standalone atomic job or extra stage is added. Native PModel atomic-density work and surviving artifact pins retain the original accounting; total internal solver calls are not inferred from file count. All method/state/gradient admissions remain unchanged and successful core status remains pending method review.

## Interpretation and remaining gates

The change directly addresses16611's documented3100MB-versus4650.8MB per-rank shortfall and OOM without changing the requested Hamiltonian. Fewer ranks reduce replicated memory and alter the runtime allowance; they can also alter floating-point reductions and the SCF trajectory. Neither memory sufficiency nor convergence is proven by this static review, and a6200MB MaxCore is not a hard process-RSS limit.

Actual zero-call preparation must pass before release. Actual grant/runtime, PModel/rotation, SCNL, full native method/state/gradient evidence and costs still need review after execution. Preserve16602 and16611 unchanged. This one declared recovery provides no automatic resource escalation, SCF fallback, displaced pair, optimization, physical-state or force-accuracy qualification.
