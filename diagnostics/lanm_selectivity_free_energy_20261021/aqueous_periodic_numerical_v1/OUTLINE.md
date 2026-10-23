# Finite aqueous periodic numerical qualification

**Source-only proposal for root integration; no molecular call, submission or MD release.** Test the actual11958atom aqueous systems and the frozen17605 water model in the unchanged BIO water/salt environment. This is numerical implementation qualification, not evidence that the Na/Cl interactions, hydration thermodynamics or ensembles are physically accurate.

**Execution held.** Root reports that the actual nine-water transfer job17670 fails three of four checks: paired energy error −4.951 kcal/mol and both individual force checks fail. This outline remains an unbound design. The frozen candidate is not ready for bulk sampling; no periodic inputs, native calls or MD are authorized by this document. The matrix below is prospective and requires a subsequent root decision after the coordination-transfer diagnosis.

## Exact starting authorities

- Mapping manifest: `workspaces/lanm_selectivity_free_energy_20261021/aqueous_water_mapping_v1/prepared_v1/MANIFEST.json`, SHA`0018cc6015898d003978bf287029e1f3532c39ae2edba0f95f18671de1220ea6`. Use only its two `aqueous_endpoints`:11958atoms,3980waters,Na7/Cl10,netzero,49.324Å cubic cell, metal native index1/type900or901, common mass138.91. Nonmetal coordinates, types, bonds and BIO bytes stay exact. The independent mapping review is `aqueous_water_numerical_v1/independent_review_v1/REVIEW.json` SHA`ad27ce03`.
- Require actual admission of the separate two-cell SI→BIO aquo8 numerical mapping check before periodic calls. Its analytic Urey–Bradley accounting must preserve BIO1.5537Å; do not change it to SI1.5337Å to obtain agreement. No repeat of those two cells is part of this inventory.
- Reuse the already qualified double-precision CUDA12/H200 engine, `gpu_engine_feasibility_v1/build_v6_cuda12/tinker9`, SHA`18cbeaf895f8bfc25a833e7e91354111ba5321f10e80b1aeb0a1d087a3f8523f`. Retain its pinned CUDA12 runtime libraries and build manifest. CPU analytic engine: `explicit_model_v1/tinker_build/testgrad`, SHA`d82b2872136785cadf0c0763aafef60ba1bdd89283f685850485a07a875f6dc2`.
- Existing GPU parity result `gpu_engine_qualification_v1/execution_v2/RESULT.json`, SHA`0e0b3acde13b537543aa87ef5a342c71fe3b0207dc955908e88250b779bf497a`, passed four POLPAIR Ewald/non-Ewald native cases (182 assertions) and four small real fixtures. Existing vacuum path result `alchemical_path_v1/run_v1/RESULT.json`, SHA`a5f7976e326e11248eb2d6408610c6f2cacbe95a3a02a2390f8cf138eb546ca4`, passed20 evaluations. Reuse those records; do not repeat their suites. Neither proves this periodic box correct.

## Periodic key and actual readback

Root should freeze the following finite numerical policy before implementation/run, consistently in CPU and GPU inputs. It follows the existing paired-leg pilot design while retaining the mapped source polarization tolerance:

```text
ewald
ewald-cutoff 9.0
ewald-alpha 0.4
pewald-alpha 0.4
pme-grid 64 64 64
pme-order 5
ppme-order 5
vdw-cutoff 9.0
vdw-taper 8.1
polarization MUTUAL
polar-eps 1.0e-8
polar-iter 100
digits 10
openmp-threads 1
verbose
```

Replace the existing identical scalar controls in the mapped key rather than accumulate conflicting copies. Preserve all model records. Add `polar-print` for CPU and `debug` for GPU so the actual final induced-dipole residual is emitted. These are logging controls, not a changed Hamiltonian. Do not set `EWALD-BOUNDARY`: both pinned native sources default to TINFOIL, and the presence of that keyword selects VACUUM regardless of a following word. Omit `VDW-CORRECTION`; the finite policy has no dispersion tail correction. No ALPB, restraints, constraints, frozen atoms, hydrogen-mass repartitioning or dielectric substitution.

Read back atom/type/class/charge inventory, box dimensions, active bonded/permanent/polarization/VDW terms, endpoint pair assignments, mutual method, polarization tolerance/cap, Ewald flag, permanent/polarization alpha0.4,64³ grid and both order5 values. `kewald.f` prints the PME table under verbose. Record actual CPU FFT implementation and GPU double-precision CUDA device/runtime. Pin source-resolved tinfoil/default-no-tail and keyword semantics alongside the key; do not falsely label an unprinted setting as observed native readback. If a material runtime field cannot be observed or unambiguously established from the exact parsed input and pinned source, leave that field unavailable and stop qualification for review; no energy rerun merely to hide missing evidence.

Important source detail: CUDA AMOEBA `src/cu/amoeba/pcg.cu` prints its final residual only when `debug` is enabled. It stops at `politer` and has a separate hard100-iteration failure branch. Therefore require the printed final RMS residual below the requested `polar-eps` for every call, with finite values and completed output; zero exit and absence of a warning are insufficient. Use the same positive residual admission on CPU. At a printed boundary, preserve precision uncertainty; do not infer a pass from rounded output.

## One proposed finite inventory:24 native programs maximum

Root declares the matrix and thresholds before execution. Execute stages in this order; failure stops subsequent stages and retains the full denominator. Each program performs one analytic E/F calculation. There are **no native Cartesian finite-difference requests**.

| Stage | Inputs | CPU programs | GPU programs | Purpose |
|---|---|---:|---:|---|
| A | Exact physical La900 and Dy901 boxes |2|2|Periodic endpoint E/F parity and positive polarization/readback admission |
| B | One fixed matched box, artificial path902 at λ=0,.25,.5,.75,1,.4999,.5001,.49995,.50005 |1 at.5|9|Endpoint identity, middle CPU/GPU parity, finite continuous parameter path and two-step λ derivative |
| C | Same midpoint on both engines, only polar-eps1e−9 |1|1|Declared tighter-polarization sensitivity; original1e−8 rows reused |
| D | Midpoint at two fixed spatial directions, each q=±1e−4 and±5e−5Å |0|8|Two directional periodic force/energy consistency checks |
| **Total** | |**4**|**20**|**24 programs; no fit or dynamics**|

StageD direction1 moves only metal index1 along+x. Direction2 rigidly translates the first complete water (native indices2,3,4) along+x, retaining its O–H/H–H geometry. Freeze those indices and directions from the prepared topology before any energy is read. For either direction, compare the central energy secants to the analytic derivative `−Σ F_i,x` over exactly the moved atoms. Two step sizes test local truncation/roundoff sensitivity; they do not check every Cartesian component. Do not request `testgrad ... Y Y` on the full box: that would trigger71,748 internal displaced energies per evaluation.

No new CPU/GPU `analyze` runs are required by this24-program proposal because analytic `testgrad` already returns total E/F and the verbose initialization/readback. If root needs native energy components beyond that output, declare their extra program count separately before launch; do not silently add them to the cap.

## Minimal mapped parameter path

Keep the existing source `alchemical_path_v1/PATH.py` unchanged. A thin new renderer must check that type/class902 is unused in the complete BIO plus900/901 namespace, then implement the same defined interpolation in that namespace. Use fixed charge+3, zero higher metal multipoles, mass138.91 and no metal bonds/groups. Atomic number57 is bookkeeping, not a fractional quantum element.

Resolve both mapped endpoints first for every represented nonmetal type402,403,405,414 (classes73,74,76,85). Interpolate metal polarizability and each resolved pair Rmin,epsilon and Thole linearly. Explicitly write four VDWPR class pairs and four POLPAIR type pairs with902, preserving all native mixing at λ0/1, including water-H/Na/Cl. Interpolate the intrinsic metal values for its isolated definition as in the prior path; do not drop or scale a pair override toward zero. Every evaluation solves mutual induced dipoles afresh. The exact physical endpoint rows from stageA are the authority for artificial λ0/1; no old vacuum energy can substitute.

## Commands and admission

Run each call in a fresh directory containing pinned `system.xyz`, `system.key`, `amoebabio09.prm`. Fully resolved executable paths are those pinned above:

```text
CPU: testgrad system.xyz Y N
GPU: tinker9 testgrad system.xyz Y N
```

Store stdout/stderr, argv, return code, allocation/device/runtime, input/output hashes, wall and CPU/GPU costs. Use unchanged `alchemical_path_v1/RUN.py::parse_analytic` (source`fdd51968`) with the existing failure/finite-number helper from `gpu_engine_qualification_v1/RUN.py` (`991fe5e3`). Require exactly one total energy,11958 ordered analytic atom rows/35874 finite components, no Numer rows, and the completed total-gradient-norm summary. Add the periodic/positive-residual readback above; do not loosen the old parser to accept truncated output. Preserve producer and request provenance separately.

Retain previous numerical tolerances prospectively: CPU/GPU and artificial endpoint equality `|ΔE|≤1e−4kcal/mol`, max force-component difference≤1e−3kcal/mol/Å. StageC must meet those same E/F limits versus each engine's own1e−8 midpoint, and pass cross-engine comparison. Both StageD secant-to-analytic discrepancies and their two-step difference≤1e−3kcal/mol/Å. The two midpoint λ energy secants must differ≤1e−3kcal/mol; these are finite-difference derivatives, not analytic `dU/dλ`. Require every scalar and every separate La/Dy endpoint gate; paired cancellation cannot rescue an endpoint. Preserve any precision-limited result as unresolved, with no automatic step/tolerance change.

One root-owned standard_expand H200 allocation with one GPU and four host CPUs can run the GPU calls and the four one-thread CPU controls sequentially, reusing the successful runtime arrangement. All numerical thread controls remain1; default scheduler memory/time. Record actual costs before any larger campaign. No login GPU, source mutation, automatic native fallback, retry, minimization or trajectory is included.

## Source evidence and limits

Inspected the pinned GPU tree's `tinker/source/kewald.f` (`fad6623e`; defaults/keyword parsing lines71,150–184; verbose PME table279–290), `cutoffs.f` (`b7625c69`; Ewald/cutoff/taper semantics), `prmkey.f` (`56207706`; POLAR-ITER/EPS and VDW-CORRECTION), and CUDA `src/cu/amoeba/pcg.cu` (`80a3467a`; residual/stop/report159–183). Corresponding CPU `tinker_source/source/kewald.f` has the same key/boundary behavior. `control.f` recognizes DEBUG/VERBOSE. The mapping and prior build authorities retain full hashes; this outline changes no source.

Passing would qualify these numerical routes and this static periodic parameter path at two starting endpoints and one middle state. It would not validate salt association, bulk hydration free energy, coordination distributions, time integration, equilibration, box-size/cutoff convergence, λ-overlap or exchange. Explicit neutrality does not remove all finite-size effects. The physical aquo9 screen is separate and is owned by Planck; its current failure is recorded above. Neither numerical parity nor a single stable box substitutes for physical transfer. Protein donor validation and any MD release remain separate.
