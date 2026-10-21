# Independent terminal review:16611

21 October 2026 UTC. **Confirmed memory failure during the first molecular Fock build. PModel initialization and the requested orbital rotation succeeded. No molecular energy/gradient result is admitted.** A separately declared resource-only recovery is supported; an unchanged-resource retry is not. Review used saved output/JSON, source and hashes only, without scientific replay or jobs.

## Exact evidence

| Artifact under run/execution unless noted | SHA256 |
| --- | --- |
| RESULT.json | `0fc5e6a1ddd98de6bfeb233d0f385f315e6f199a9e608e31009b9d3371119d3c` |
| central_La/endpoint.out | `9d877978f816272fe22a2d03bc876a843378e5dc997df501cb34cbedb51ed9c3` |
| central_La/endpoint.out.execution.json | `1a073ff0296d8277a00c950d8c7bd3f656481a824c58002767647c3149bcad2e` |
| NATIVE_EVIDENCE.json | `fc62b81b9babed6a7de9e45d0abd63fdffec9cfaf4f2f157ad2811f6ddf3ff23` |
| ATOMIC_GUESS_ARTIFACTS.json | `f5eafdb317feb5d32ba37f3426294adc6f0fefc452e2e564ceea5ddb3bdc0f4d` |
| central_La/endpoint.runtime_atom57.out | `9628688e395b22c55c1ff1a069a8c72f7a45dbdb05bc8093373ca0cd09ebd5af` |
| run/slurm-16611.err | `e305a4e3ee4e5275352ada997516b0cd3c800b23df0cfc68514b12500a5527a1` |

Selected source/configuration/manifest, all non-null runner artifacts, external basis and all five surviving atomic-guess artifact hashes match. Source remains `27ac4e80…`, run configuration `b654d846…`, input `4359537c…`, XYZ `55e07d81…`, basis `7389c56f…`. The result remains `invalid`, comparison unavailable, with one runner/whole-molecule attempt and zero other requested endpoints or engines. Native receipt normal termination, SCF convergence and expected engrad are false. As in16602, top-level ORCA returncode0 does not override its internal failure. No engrad exists.

## Why this is an OOM, not the earlier guess failure

Native lines2155–2163 show MODEL POTENTIAL and successful return from La/ECP atomic-density generation. Lines2181–2186 print the514(1)/515(1) rotation with the requested20-degree cosine/sine and completed initial guess. Atomic multiplicity2 is an initialization detail, not a change to the declared molecular singlet.

LEAN-SCF then warns at2204–2206 that3100.0MB is available versus4650.8MB estimated. It reaches the first incremental Fock formation at2212, without a printed molecular iteration. MPI reports rank0 killed by signal9 and ORCA error termination in LEANSCF. Slurm stderr explicitly records an `oom_kill` in16611.batch. The terminal event records OUT_OF_MEMORY, exit0:125,962s and24 allocatedCPUs. Thus missing terminal SCF residuals are expected, not evidence of failed convergence or a parser-only problem.

The actual grant was24MPI/24CPUs/96000MiB with one thread per rank and MaxCore3100MB. The earlier startup matrix-fit check and low Python-parent RSS did not characterize peak molecular SCF memory. Initial-density or temporary energy-file values are not a converged molecular result. SCNL/libXC, external La ECP46, UKS state and the declared scientific settings remain as observed previously; no physical-state or force qualification follows.

## Resource-only recovery assessment

The pinned renderer `46356652…` distinguishes allocated CPUs from local MPI slots, limits ranks by both and derives per-rank MaxCore from granted memory divided by workers×ranks. The pinned runner `b9278e74…` and base renderer allow fewer ranks than allocated CPUs while retaining explicit one-thread settings. The current TARGET allocation guard, however, requires both granted CPUs and slots to equal24. Any different resource layout needs a separately versioned exact allocation/rank guard and corresponding release, not an environment override of the frozen adapter.

Root initially considered24MPI with48 allocatedCPUs/default192000MiB. That arrangement is supported by the generic renderer in principle, but requires scheduler permission beyond the earlier observed24CPU/job QoS limit; it must not be assumed legal.

The subsequently reported bounded recovery design instead uses **12MPI slots ×2 allocatedCPUs/task =24 allocatedCPUs, default96000MiB, unchanged25% reserve and derived MaxCore6200MB**, retaining one numerical thread per rank. This resource layout is coherent: it raises the per-rank allowance above the printed estimate and reduces replicated rank memory while keeping the existing CPU allocation ceiling. These resource numbers are the declared proposal, not an executed recovery or a guarantee against later peak memory.

Preserve the exact molecular input/XYZ/basis bytes, PModel and beta rotation, SCNL behavior, all method/state/grid/SCF controls and finite native admission. Change only explicitly declared allocation/MPI/runtime-memory handling. Keep both16602 and16611 failures immutable, use fresh unseeded files and no checkpoint/GBW fallback. Actual resource observation, native outcome and total cost must govern the single recovery; no automatic escalation, solver change, displaced pair or optimizer follows.
