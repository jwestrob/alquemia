# Installed GPU portability: exact AMOEBA EXCHPOL route is absent

**Direct use on the installed Tinker9 GPU binary is blocked.** S2U exchange-polarization kernels exist, including periodic geometry and gradient/virial code, but they are dispatched only through HIPPO polarization. The CPU diagnostic uses unchanged AMOEBA/Thole parameters with no charge penetration. Source inspection predicts that its EXCHPOL response would be bypassed on the current GPU route; this has not been measured by a native call.

| Capability | Evidence-based status |
|---|---|
| EXCHPOL parameter loading and S2U overlap | Present. `src/hippo/expol.cpp` copies class-resolved arrays and recognizes S2U. `include/seq/pair_alterpol.h:16` uses the same t=0.5√(bᵢbⱼ)r, squared overlap and radial derivative as the CPU. |
| Tensor inversion and exchange derivative CUDA kernels | Present in source and linked binary. Read-only `nm -C` finds `alterpol_cu`, `dexpol_cu`, and `induceMutualPcg4_cu` with double-precision arguments. |
| Current AMOEBA energy/force path | **Absent dispatch.** [energy.cpp](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/gpu_engine_feasibility_v1/tinker9_source/src/energy.cpp:158) rejects HIPPO polarization unless `use_chgpen` is true. [AMOEBA induce](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/gpu_engine_feasibility_v1/tinker9_source/src/amoeba/induce.cpp:108) always uses two-response PCG1, without EXCHPOL tensors. Its non-Ewald and Ewald polarization paths both call this function; the fused CUDA `emplar` route does too. No `dexpol` call is present in those AMOEBA paths. |
| HIPPO exchange induction/derivatives | Present. [induce2](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/gpu_engine_feasibility_v1/tinker9_source/src/hippo/induce.cpp:40) dispatches tensor-aware PCG4; [HIPPO epolar](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/gpu_engine_feasibility_v1/tinker9_source/src/hippo/epolar.cpp:235) adds `dexpol`. PCG4 calls **charge-penetration** permanent/mutual fields, so routing our keys there is not a faithful shortcut. |
| Periodic support for the exact declared Hamiltonian | **Unqualified and blocked by the same dispatch gap.** HIPPO CUDA kernels accept image/box arguments and a repulsion switch; tests include periodic salt/water examples. This is source capability, not an executed test of our AMOEBA term. |

The actual prepared q=.05 key and parameter file contain no active `chgpen` record. Adding one to reach the HIPPO route would change the Hamiltonian. Parsing EXCHPOL, allocating its arrays, or finding kernel symbols does not prove that those arrays participate in the selected energy/gradient path.

## Installed artifact identity

- Source HEAD: `44bcd7c898e68f6827c112df79b8361e443730d2`. All decisive physics files inspected here match HEAD. The current overall working diff only marks the pre-existing Catch2 submodule change; no engine files were edited in this review.
- GPU binary: `/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/gpu_engine_feasibility_v1/build_v6_cuda12/tinker9`, SHA `18cbeaf895f8bfc25a833e7e91354111ba5321f10e80b1aeb0a1d087a3f8523f`; matches the bound qualification config. Build cache: CUDA, HOST=OFF, PREC=double. Existing H200/POLPAIR qualification does not establish this missing AMOEBA dispatch.
- Exact CPU probe manifest: `2f22ff157b09118c8f36e6abccadb4ecaf6953c1feecef1719b75a7937ec8f2a`.

Decisive inspected file SHA256 pins:

| Source relative to the Tinker9 tree | SHA256 |
|---|---|
| `src/energy.cpp` | `27a8e630fd741aeba2b9fac7a385fe77e00bf2add40ba0da7c4f206d1ee82297` |
| `src/amoeba/induce.cpp` | `63a1370dda5adfc1561a36de582bab5d5c7d309b22cebe00e2d813c63008e020` |
| `src/amoeba/epolar.cpp` | `5d9fce71954fd03580d9046113c804f15cf6ea7198ddfb860d6fe81b9e37ed6f` |
| `src/hippo/induce.cpp` | `7a2888f6e279e56276a105cdf341f128f0d302334686f5d0688dcdaff87adf94` |
| `src/hippo/epolar.cpp` | `150e312f3bc33102bd0ea31932f58e2b808d95966930ef021ad63ba63e0b3f39` |
| `src/hippo/expol.cpp` | `59ad7941c89ea8c89205024e98cecf20aa9f4d8eb2032fe439a35e5ad382b1cb` |
| `src/cu/hippo/expol.cu` | `0be1ed51308b68a1a01e57deb75e4268586a127be88c47007758d87402b0ffed` |
| `include/seq/pair_alterpol.h` | `92fddd37fe269da4d12998d994fc8df9ce40c81b82b799d1f43a0d638d5cd63c` |

## Smallest conditional numerical gate

Only if the CPU mechanism proves useful, first preserve AMOEBA's permanent and Thole fields, two-response induction, predictor, force and virial conventions in a separately reviewed implementation. Reusable GPU tensor/S2U kernels reduce that work, but this is more than a keyword or branch toggle. No such edit/build is authorized or performed here.

1. **Vacuum route gate: six GPU program calls.** Reuse CPU18295 outputs and exact prepared inputs for q0/q.05 W9 La/Dy (four cells), plus the metal-absent W9 shell at q0/q.05 (two). Run `tinker9 testgrad system.xyz Y N`; for the two q.05 metal clusters replace `N` with `Y 0.00001`, giving exactly **336 internal displaced-energy evaluations**. Compare each E/F and q-dependent change with CPU, using inherited |ΔE|≤1e−4 kcal/mol and maximum |ΔF|≤1e−3 kcal/mol/Å; require nonzero CPU response to be reproduced, zero/pure-water invariance and admitted induction residuals. Retain off-diagonal response/FD evidence. No new CPU calls are needed if those six CPU rows are admitted.
2. **Periodic static gate only after that passes:** use both already-prepared explicit-neutral aqueous endpoints, unchanged coordinates/BIO water and explicit S2U class mapping, q0/q.05 each. Map the exact diagnostic 17676 metal-water pair values; do not substitute the older 17605 overrides currently carried by the prepared boxes. Four inputs × CPU/GPU = **eight analytic-only E/F program calls**, matched Ewald, PME grid, cutoff/taper, boundary convention and readback. Verify every site's alpha before invoking this CPU derivative route. Compare total and q-change E/F under predeclared tolerances; retain full finite/induction checks. **Never request all-Cartesian finite differences on the 11,958-atom boxes.** This is a prospective inventory, not a bound panel or release. Virial/energy-conservation validation remains necessary before MD, and salt/physical transfer remains unqualified.

No molecular programs, GPU tests, builds, submissions or engine edits occurred. One local thread; only source, saved metadata and binary symbol/hash reads. The CPU mechanics diagnostic remains interpretable independently of this deployment limitation.
