# One memory-only recovery: independent review

**Conclusion:** one explicitly recorded recovery at 256 GiB is scientifically consistent with the existing La-aquo compatibility pilot, provided the deck, executable, assets, 12 MPI × 2 OpenMP layout, algorithms, thresholds and admission gates remain exact. The evidence identifies an allocation failure, not a demonstrated ALMO/PBE0/ECP incompatibility. The proposed memory is reasonable headroom, not a proven sufficient peak-memory bound. This review authorizes no submission and changes no source or input.

## Actual failed attempt

Job 21618 is terminal OUT_OF_MEMORY, with 126 scheduler seconds and 24 allocated CPUs in the saved completion receipt. RUN_RESULT records native exit 9 after 118.65459037292749 seconds; the whole worker took 126.61188986198977 seconds. The launcher records signal-9 termination. Root reports MaxRSS 98,304,464 K against 96,000 MiB; that specific MaxRSS value was supplied by root, whereas the terminal state, elapsed values and raw output were independently read here. It should remain joined to root's actual accounting evidence, not interpreted as a per-MPI-rank measurement.

The raw output contains exactly one `MOLECULAR GUESS:` banner, for fragment 1 (La, charge +3, multiplicity 1). There are zero fragment `SCF run converged` markers, zero `Single-molecule energy:`, zero `Energy of the initial guess:`, zero ALMO DIIS/PCG rows, zero final `ENERGY| Total FORCE_EVAL`, and no normal program-end record. The output does contain six auxiliary **neutral-La atomic-guess** iterations and their energy. Those are initialization data, not the requested isolated M3+ fragment energy, R/F/P/U, or a full-system electronic result. No scientific cell has completed and none can be reused as an admitted ALMO result.

The original unsuccessful execution and parser failure must remain unchanged. A new isolated execution directory permits a fresh, unchanged initial guess without confusing partial files with an admitted restart. Count the original molecular-parent attempt and the recovery separately; a recovery success would not erase the first failure or its cost.

## Why 256 GiB is a defensible resource choice

The actual output prints two sets of large grids: the parent system and the currently initialized fragment. Each includes a 900 × 900 × 900 MT auxiliary grid and a 625 × 625 × 625 charge grid, distributed over 12 processors. These sizes are consequences of the unchanged 30-A cell and tight declared cutoff, not a new electronic method. The largest grid has 729 million points: a single real-double field alone is about 5.43 GiB globally, or about 10.86 GiB if represented by one complex-double field. These are illustrative storage units, **not** a sum of simultaneously live allocations or a peak-memory prediction.

The source supports the observed memory pressure:

- `src/mscfg_methods.F:239–306,348–349` keeps the parent environment while creating a separate fragment environment with the same full cell; that fragment is released only after its energy and orbitals are obtained. Small atom count therefore does not imply a small molecular-guess grid.
- `src/pw/mt_util.F:72–103` creates a higher-cutoff auxiliary pool and real/reciprocal temporary fields for the nonperiodic MT screen. `src/pw/pw_types.F:117–133` allocates arrays using actual local grid sizes.
- The raw real-space output shows substantial halos on the 625-point grid (border 74, local plane maxima 357 and 305). Several smaller real-space grids are fully replicated. Thus one cannot estimate total memory from the largest distributed field alone or divide all storage simply by 12.

256 GiB is roughly 2.73 times the original 96,000-MiB allocation. This is meaningful headroom over the observed cgroup failure boundary. It is not an empirically established ceiling: the shell fragment, full SCF, exchange work, FFT buffers and later temporaries had not yet completed. This source review cannot promise either memory sufficiency or completion within the remaining time. No lower cutoff, smaller cell, altered MPI layout, screened/truncated exchange, modified basis, or method substitution is implied.

## Scope and fixed budget

Retain the exact input SHA `15b07a8ba0a5cb8037847908dd87c6977fa50ed53c8a07ca7489fc7a5d92883b`, the existing executable/assets and all settings. Retain the declared 1e-4-Eh cross-code full-energy gate and 1e-8-Eh arithmetic/order gates; actual fragment and ALMO convergence remain independently required. The other three quartet cases remain gated by the completed pilot.

Do not reset the original 7,200-second native budget. Deducting the recorded first native elapsed leaves at most 7,081.345409627073 seconds (an integer cap of 7,081 seconds). A more conservative deduction of the full first worker elapsed leaves 7,073.388110138010 seconds (integer cap 7,073 seconds). Root should state which accounting convention is used; neither increases the original budget. The saved scheduler allocation is 126 × 24 = 3,024 CPU-seconds, distinct from RUN_RESULT's worker-wall estimate of 3,038.6853566877544 allocated CPU-seconds.

Current AGENTS.md explicitly routes subsequent work through standard_expand. Its no-memory-override wording belongs to the labeled historical OMol16680 section; it does not impose a new blanket restriction on this separately declared CP2K pilot. Root owns resource availability, binding and submission. A failed recovery must remain an explicit failure; this review supplies no automatic further attempt.

## Evidence pins and review scope

All paths below are under `workspaces/lady_electronic_decomposition_20261023/cp2k_quartet_v1/tasks/aquo8_central_La/`:

| Evidence | SHA256 |
|---|---|
| `RUN_RESULT.json` | `e2da1f0beef956cc53b8744982b5ee1db3d99f574c61a9db404d14ee7392ce1c` |
| `RESULT.json` | `b1ab02b1ddceb66ebd163f1f8342d696b54921db38862b5c8588fe0af6148f8b` |
| `completion_event.json` as read | `225dbf0bd402585689c6483a8831a0abe52d10679f4c5df03ef1a08002827ed4` |
| `execution/endpoint.out` | `d997a2053481a39285fa59e10fc93c558a839b002bbd01468bec2dec63fcf23c` |
| `execution/native_launcher.log` | `e13b4d1e86f135b7f63a04475d0f01290872233aca4df298a946dc992f37da98` |

Source authority remains the pinned actual CP2K 2026.2 release/archive described in `cp2k_inputs_v1/PILOT_SETTINGS_REVIEW.md`, including the already reviewed Dy reader patch. This review used saved text, source reads, hash/count metadata and elementary storage/budget bookkeeping. It performed no scientific calculation, runtime import, build, input generation or scheduler action.
